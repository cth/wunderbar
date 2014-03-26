require 'progressbar'

class Matcher
	# Initialize with 
	def initialize(config) 
		@config = config
		# Read all barcoding genotypes from plink files
		@barcode_genotypes = Genotypes.new(@config["barcode-file"])
		# Read only barcoding genotypes from testset plink files
		@test_genotypes = Genotypes.new(@config["file"],@barcode_genotypes.snps)
		puts "barcodes: " + @barcode_genotypes.to_s
		puts "test snps: " + @test_genotypes.to_s
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
			#puts "snp=#{snp} matching=#{matching} mismatch=#{mismatch} unknown=#{unknown}"
		end
		pbar.finish
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
				if pvalue < @config["match-p"].to_f and evalue < @config["match-e"].to_f
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

				next unless @mismatches.include?(indv_bar) 

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
				next if p_value > @config["match-p"].to_f
					
				inf_factor = ld_inflation_factor(mismatch,matching,p_value) unless @config["mcmc"].nil?
				p_value = p_value * inf_factor 

				if p_value <= 0.05 then
					@swaps << { :original_label => indv_test, :new_label => indv_bar, :pvalue => p_value , :ld_inflation_factor => inf_factor , :evalue => mismatch.size*p_value, :match => matching.size, :mismatch => mismatch.size, :unknown => unknown } 
				end
			end
		end
		pbar.finish
	end

	def create_swaps_report(out_file)
		if @swaps.size == 0 then 
			puts "No swaps detected"
		else 
			File.open(out_file, "w") do |file|
				keys = @swaps.first.keys
				file.puts (keys.map { |x| x.upcase }).join("\t")
				swaps = @swaps.sort { |a,b| a[:pvalue] <=> b[:pvalue] } 
				swaps.each { |swp| file.puts (keys.map { |k| swp[k] }).join("\t") }
			end
		end
	end

	# Some SNPs may in close LD, in which case we would expect
	# an overrepresentation of certain patterns of genotypes in @test_genotypes 
	# p.values are calculated assuming independence. This fn calculates a factor to adjust for this.
	def ld_inflation_factor(mismatches,matches,p_value)
		avg_inflation_ratio = 1
		if mismatches.size > 0 then
			snp_subset_sizes = []
			# Create sample distribution of binomials
			mismatches.size.downto(2) { |x| snp_subset_sizes << x }
			inflation_ratios = []
			# Scale number of trials relative to extremity of p-value
			# FIXME: 100 should be a configurable factor
			1.upto(@config["mcmc"].to_i) do |i|
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
			avg_inflation_ratio = inflation_ratios.inject { |a,b| a+b } / inflation_ratios.size
			puts avg_inflation_ratio
		end
		avg_inflation_ratio
	end
end
