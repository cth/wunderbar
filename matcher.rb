require 'progressbar'

class Matcher
	# Initialize with 
	def initialize(config) 
		@config = config
		# Read all barcoding genotypes from plink files
		@barcode_genotypes = Genotypes.new(@config[:barcode_file])
		# Read only barcoding genotypes from testset plink files
		@test_genotypes = Genotypes.new(@config[:file],@barcode_genotypes.snps)
		puts "barcodes: " + @barcode_genotypes.to_s
		puts "test snps: " + @test_genotypes.to_s

		puts "Sampling p-value distributions for mismatch binomials:"
		sample_mismatch_distributions( @config[:mismatch_adjust].to_i )
	end

	# returns true if value is nil or zero
	def nilz(x)
		x.nil? or x==0 
	end

	# Match each SNP to its equivalent SNP in other dataset
	# This is to check that snp-wise correspondance is ok
	# If this is ok, we can then use to detect outlying individuals
	def match_snps
		@snp_stats = {}
		@mismatches = {}

		pbar = ProgressBar.new("Single barcode match", @barcode_genotypes.snps.size)
		@barcode_genotypes.snps.each do |snp|
			pbar.inc
			matching,mismatch,unknown = 0,0,0
			@test_genotypes.each_by_snp(snp) do |indv,gt|
				if not @barcode_genotypes.include_individual?(indv) or nilz(gt) or nilz(@barcode_genotypes[indv,snp]) then
					unknown = unknown + 1
				elsif gt == @barcode_genotypes[indv,snp] 
					matching = matching + 1
				else
					mismatch = mismatch + 1
					@mismatches[indv] ||= []
					@mismatches[indv] << snp
				end
			end
			@snp_stats[snp] = { :match => matching, :mismatch => mismatch, :missing => unknown, :total => matching+mismatch+unknown }
		end
		pbar.finish
		# calculate e-values and p-values
		@individual_match = {}
		@test_genotypes.individuals.each do |indv|
			pvalue = 1
			if @mismatches[indv].nil? then
				evalue = "NA" 
			else
				@mismatches[indv].each do |snp| 
					pvalue = pvalue  * (@snp_stats[snp][:mismatch].to_f / @snp_stats[snp][:total].to_f)
				end
				# Ok, this is messy. We want to take into account that some SNPs may be missing.
				avg_samples = (@mismatches[indv].map { |snp| @snp_stats[snp][:total] }).inject { |r,v| r+v } / @mismatches[indv].size.to_f
				evalue = pvalue * avg_samples
			end
			@individual_match[indv] = { :pvalue => pvalue, :evalue => evalue } 
		end
		puts @individual_match.inspect
	end

	def create_match_report(output_file)
		File.open(output_file,"w") do |file|
			file.puts "ID\tNMISMATCH\tERATE\tP-VALUE\tE-VALUE"
			(@mismatches.to_a.map { |k,v| [v.size,v,k] }).sort.reverse.each do |numsnps,snps,particid| 
				pvalue = 1
				# Calculate a p-value for particid (probability of seing this particular combination of mismatches at random) 
				snps.each { |snp| pvalue = pvalue * (@snp_stats[snp][:mismatch].to_f / @snp_stats[snp][:total].to_f) }
				# In the same breath, calculate average number of samples for the compared snps 
				avg_samples = (snps.map { |snp| @snp_stats[snp][:total] }).inject { |r,v| r+v } / snps.size.to_f
				# Calculate e-value:
				evalue = avg_samples * pvalue
				if pvalue < @config[:p_max] and evalue < @config[:e_max]
					file.puts "#{particid}\t#{numsnps}\t#{numsnps.to_f/@barcode_genotypes.snps.size}\t#{pvalue}\t#{evalue}"
				end
			end
		end
		puts "#{@mismatches.size} samples with barcoding errors written to #{output_file}"
	end

	def swap_match
		# Construct alignment matrix between all mismatching samples
		@swaps = []
		pbar = ProgressBar.new("Swap matching", @barcode_genotypes.individuals.size * @test_genotypes.individuals.size)
		@barcode_genotypes.individuals.each do |indv_bar|
			@test_genotypes.individuals.each do |indv_test|
				pbar.inc
				next if indv_bar == indv_test

				#next unless @mismatches.include?(indv_bar) 

				unknown = 0
				matching = []
				mismatch = []
				@test_genotypes.each_by_indv(indv_bar) do |snp,gt|
					other_gt = @barcode_genotypes[indv_test,snp] 

					if nilz(gt) or nilz(other_gt) then
						unknown = unknown + 1
					elsif gt == other_gt then
						matching << [ snp, gt ]
					else
						mismatch << [ snp, gt ]
					end
				end

				total = matching.size + unknown + mismatch.size 
				# p-value for match: Product of all the genotype frequencies (in @test_genotypes) for matching SNPs 
				p_value = 1
				matching.each { |snp,gt| p_value = p_value * @test_genotypes.snp_freq[snp][gt] }
				#mismatch.each { |snp,gt| p_value = p_value * binomial(@barcode_genotypes.snps.size,mismatch.size) }
					
				p_value = p_value * Binomial.new(matching.size + mismatch.size, matching.size).to_i

				inf_ratio = ld_inflation_ratio(mismatch,matching,p_value)
				p_value = p_value * inf_ratio
				
				p_value = (p_value > 1) ? 1 : p_value # cap p-value at 1
 
				e_value = mismatch.size*p_value

				posterior_prob = (1-p_value) * (@config[:peturb]+(1-@individual_match[indv_bar][:pvalue])) * (@config[:peturb]+(1-@individual_match[indv_test][:pvalue]))

				next if (p_value > @config[:p_max].to_f or e_value > @config[:e_max])

				@swaps << { :original_label => indv_test, :new_label => indv_bar, :posterior => posterior_prob, :pvalue => p_value , :ld_inflation_ratio => inf_ratio , :evalue => e_value, :match => matching.size, :mismatch => mismatch.size, :unknown => unknown } 
			end
		end
		pbar.finish
	end

	def create_swaps_report(out_file)
		if @swaps.size == 0 then 
			puts "No swaps detected"
		else 
			keys = @swaps.first.keys
			swaps = @swaps.sort do |a,b| 
				if a[:posterior] == b[:posterior] then
					b[:pvalue] <=> a[:pvalue]
				else
					a[:posterior] <=> b[:posterior]
				end
			end.reverse
			File.open(out_file, "w") do |file|
				file.puts (keys.map { |x| x.upcase }).join("\t")
				swaps.each { |swp| file.puts (keys.map { |k| swp[k] }).join("\t") }
			end

			if not @config[:dot_file].nil? then 
				File.open(@config[:dot_file], "w") do |dot_file|
					dot_file.puts "digraph G {"
						swaps.each do |swp|
							if swp[:posterior] > 0.95 then
								dot_file.puts '"' + swp[:original_label] +  '" -> "' + swp[:new_label] + '"'
							end
						end
					dot_file.puts "}"
				end
			end
		end
	end

	# Some SNPs may in close LD, in which case we would expect
	# an overrepresentation of certain patterns of genotypes in @test_genotypes 
	# p.values are calculated assuming independence. This fn calculates a factor to adjust for this.
	def ld_inflation_ratio(mismatches,matches,p_value)
		avg_inflation_ratio = 1
		if mismatches.size > 0 then
			snp_subset_sizes = []
			# Create sample distribution of binomials
			mismatches.size.downto(2) { |x| snp_subset_sizes << x }
			inflation_ratios = []
			# Scale number of trials relative to extremity of p-value
			1.upto(@config[:ld_adjust].to_i) do |i|
				# randomly select number of snps in subset
				num_snps = snp_subset_sizes.random_select() 

				# Sample a corresponding number of the matching SNPs
				sampled_snps = matches.sample_n(num_snps)

				# Count number of matches for this subset 
				num_match = 0 
				@test_genotypes.individuals.each do |indv|
					# map snps to match truth value and logical and 
					if (sampled_snps.map { |snp,gt| gt == @test_genotypes[indv,snp] }).inject { |a,b| a and b } then
						num_match = num_match + 1 
					end
				end
				# We do not correct if num_match is zero
				next if num_match == 0

				# calculate probability of seeing by chance
				snp_probs = sampled_snps.map { |snp,gt| @test_genotypes.snp_freq[snp][gt] }
				next if  snp_probs.include?(nil) # FIXME: this can happen if the barcode includes a genotype not included in @test_genotypes
				expected_match_prob = snp_probs.inject { |a,b| a*b }
				expected_match = expected_match_prob * @test_genotypes.individuals.size 

				inflation_ratios <<  num_match / expected_match 
			end
			avg_inflation_ratio = inflation_ratios.sum / inflation_ratios.size if @config[:ld_adjust].to_i > 0
		end
		avg_inflation_ratio
	end

	
	def sample_mismatch_distributions(iterations)
		binomials = []
		(@barcode_genotypes.snps.size-1).downto(1) { |x| binomials << Binomial.new(@barcode_genotypes.snps.size,x) }
		@pvalue_distribution=[]
		binomials.each { |b| @pvalue_distribution[b.bottom] = sample_pvalues(b,iterations) }
	end

	def sample_pvalues(binom,iterations)
		sampled_pvalues = []
 
		pbar = ProgressBar.new("#{binom.to_s}", iterations)

		# Run sampling
		1.upto(iterations) do |i|
			pbar.inc
			# sample n snps and calculate probability of seeing by chance 
			# uniform probability for selection of each snp and genotype within snp 
			# TODO: Should it be uniform with snp??
			sample_probs = @test_genotypes.snps.map { |x| 1.to_f / @test_genotypes.snps.size }
			sampled_snps =  @test_genotypes.snps.sample_n(binom.bottom,sample_probs)
			probs = sampled_snps.map { |snp,gt| @test_genotypes.snp_freq[snp][[0,1,2,3].random_select()] }
			pvalue = probs.inject { |a,b| a*b }
			sampled_pvalues << pvalue
		end
		pbar.finish
		sampled_pvalues
	end
end
