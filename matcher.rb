require 'progressbar'

class Matcher
	# Initialize with 
	def initialize(test_genotypes_plink12, bar_genotypes_plink12)
		# Read all barcoding genotypes from plink files
		@barcode_genotypes = Genotypes.new(bar_genotypes_plink12)
		# Read only barcoding genotypes from testset plink files
		@test_genotypes = Genotypes.new(test_genotypes_plink12,@barcode_genotypes.snps)
		puts "barcodes: " + @barcode_genotypes.to_s
		puts "test snps: " + @test_genotypes.to_s
	end

	# returns true if value is nil or zero
	def nil_or_zero(x)
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
				if not @barcode_genotypes.include_individual?(indv) or nil_or_zero(gt) or nil_or_zero(@barcode_genotypes[indv,snp]) then
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
				file.puts "#{particid}\t#{numsnps}\t#{numsnps.to_f/@barcode_genotypes.snps.size}\t#{pvalue}\t#{evalue}"
			end
		end
		puts "#{@mismatches.size} samples with barcoding errors written to #{output_file}"
	end

	def swap_match
		# Construct alignment matrix between all mismatching samples
		@swaps = []
		pbar = ProgressBar.new("Swap matching", @barcode_genotypes.individuals.size * @test_genotypes.individuals.size)
		@barcode_genotypes.individuals.each do |indv_a|
			@test_genotypes.individuals.each do |indv_b|
				pbar.inc
				next if indv_a == indv_b

				unknown = 0
				matching = []
				mismatch = []
				@test_genotypes.each_by_indv(indv_a) do |snp,gt|
					other_gt = @barcode_genotypes[indv_b,snp] 

					if nil_or_zero(gt) or nil_or_zero(other_gt) then
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
				next if p_value >= 1
					
				#pvalue_adjust = mcmc_adjust(mismatch,indv_a,indv_b) 
				pvalue_adjust = 0 

				if p_value <= 0.05 then
					@swaps << { :original_label => indv_b, :new_label => indv_a, :pvalue => p_value + pvalue_adjust, :adjust => pvalue_adjust, :evalue => mismatch.size*p_value, :match => matching.size, :mismatch => mismatch.size, :unknown => unknown } 
				end
			end
		end
		pbar.finish
	end

	def create_swaps_report(out_file)
		if @swaps.nil? then 
			puts "No swap detected"
		else 
			File.open(out_file, "w") do |file|
				keys = @swaps.first.keys
				file.puts (keys.map { |x| x.upcase }).join("\t")
				swaps = @swaps.sort { |a,b| a[:pvalue] <=> b[:pvalue] } 
				swaps.each { |swp| file.puts (keys.map { |k| swp[k] }).join("\t") }
			end
		end
	end

	# Monte Carlo simulation to calculate how often
	# we get a "match" p-value at least as extreme with as many matching snps
	def mcmc_adjust(mismatch,indv)
		adjust = 0
		if mismatch.size > 0 then
			binomials = []
			mismatch.size.downto(1) { |x| binomials << Binomial.new(mismatch.size + matching.size,x) }
			binomial_probs = normalize_cnts(binomials.map { |x| x.to_i } )
			smaller_probs = []
			larger_probs = []
			trials = 100*(1.to_f/p_value).to_i
			puts "trials: #{trials}"
			puts binomial_probs.inspect
			1.upto(trials) do |i|
				# Randomly select binom.bottom barcoding snps and calculate match probability	
				# record numbers of times that it is better or worse than matching p-value respective

				binom = binomials.random_select(binomial_probs)

				# uniform probability for selection of each snp 
				sample_probs = @test_genotypes.snps.map { |x| 1.to_f / @test_genotypes.snps.size } 

				# sample n snps: 
				sampled_snps =  @test_genotypes.snps.sample_n(binom.bottom,sample_probs) 
	
				# calculate probability of seeing by chance
				prob = sampled_snps.map { |snp,gt| @test_genotypes.snp_freq[snp][gt] }
				next if  prob.include?(nil) # FIXME: potential problem
				prob = prob.inject { |a,b| a*b }

				# Put probaility into "smaller" bucket or "bigger" bucket
				if prob > p_value then
					larger_probs << prob
				else
					smaller_probs << prob
				end

				# Adjust p-value
				if smaller_probs.empty? then
					adjust = 0
				else
					fraction_smaller = smaller_probs.sum.to_f / (smaller_probs.sum + larger_probs.sum)
					mean_smaller = smaller_probs.sum / smaller_probs.size
					adjust = mean_smaller * fraction_smaller * (binomials.map { |b| b.to_i }).sum
				end
			end
		end
		adjust
	end
end
