#!/usr/bin/env ruby

require_relative 'df'
require 'rinruby'
require 'yaml'
require 'pp'

# Make sure the NCBI2R library is installed in R
R.eval "library(NCBI2R)"

# Add a method to calculate binomial coefficients to the Integer class
class Integer
	def choose(k)
		pTop = (self-k+1 .. self).inject(1, &:*)
		pBottom = (2 .. k).inject(1, &:*)
		pTop / pBottom
	end
end

class Binomial
	attr_reader :top, :bottom
	def initialize(top,bottom)
		@top,@bottom = top,bottom
	end

	def to_i
		@top.choose(bottom)
	end

	def to_s
		"#{@top} choose #{@bottom}"
	end
end

class Array
	def sum
		if self.size == 0
			0
		else
			self.inject { |r,n| r+n }
		end
	end

	def random_select(probs=nil)
		threshold = rand

		throw "probs doesnt match array length" unless probs.nil? or probs.size == self.size
		cumulative_sum = 0
		0.upto(self.size-1) do |i|
			if probs.nil? then
				cumulative_sum = cumulative_sum + (1 / self.size.to_f)
			else
				cumulative_sum = cumulative_sum + probs[i]
			end
			return self[i] if cumulative_sum > threshold
		end
		self[self.size-1]	
	end

	def sample_n(n,probs=nil)
		samples = []
		throw "Cannot samples #{n} from and array of size #{self.size}" if n > self.size
		while samples.size != n
			e = self.random_select(probs)
			samples << e unless samples.include?(e)
		end
		samples
	end
end

def normalize_cnts(arr) 
	sum = arr.sum
	arr.map { |x| x.to_f / sum }
end

class GenotypeMatrix
	def initialize
	end

	def intersect
	end
end

class BarCoder
	def initialize(barcode_calls_file,plink_stem,outputdir)  
		@barcode_calls_file, @plink_stem, @outputdir = barcode_calls_file, plink_stem, outputdir
		@particid_snps = {}
		@snp_particids = {}
		@snpinfo = {}
		@output_file_mismatch = @outputdir + "/barcode_mismatches.txt"
		@output_file_swaps = @outputdir + "/swaps.txt"
		@cutoff_detect_wrong_sample = 0.05
		@cutoff_swap_samples
	end

	def load
		if File.exists?("output/offline.yaml") then  
			@particid_snps, @snp_particids = *YAML.load(File.read("output/offline.yaml"))
			@particid_snps.each do |k,v|
				@particid_snps[k] = v.delete_if { |vv| vv == nil }	
			end	
		else
			DF.new(@barcode_calls_file, { :header => true }).each do |particid,snppos,genotype|
				puts snppos	
				if snppos =~ /(\d+)_(\d+)/
					chromosome = $1, position = $2 
				elsif snppos =~ /rs.*/ and @snpinfo[snppos].nil?
					R.eval "snp <- GetSNPInfo('#{snppos}')"
					chromosome = R.pull("snp$chr").to_i
					position = R.pull("snp$chrpos").to_i
					@snpinfo[snppos] = [ chromosome, position ] 
				elsif not @snpinfo[snppos].nil?
					chromosome, position = *@snpinfo[snppos]
				else
					throw "Strange SNP #{snppos}"
				end

				chromosome = chromosome.first if chromosome.class == Array
				position = position.first if position.class == Array
				if genotype =~ /^(A|G|C|T):(A|G|C|T)$/
					a1 = $1
					a2 = $2
				end	

				snp = {
					:particid => particid,
					:chromosome => chromosome,
					:position => position,
					:a1 => a1,
					:a2 => a2,
					:snppos => snppos,
					:snpkey =>  chromosome.to_s + ":" + position.to_s
				}
				#	puts snp.inspect

				@particid_snps[particid] ||= [] 
				@particid_snps[particid] << snp unless snp.nil? 
				@snp_particids[snp[:snpkey]] ||= [] 
				@snp_particids[snp[:snpkey]] << snp unless snp.nil?
			end
		end
	end

	def save_barcoding_data
		# Write to files: 
		File.open("output/offline.yaml", "w") do |f| 
			f.write([@particid_snps, @snp_particids].to_yaml)
		end
	end

	def write_snp_keep_list
		keep = []
		DF.new(@plink_stem + ".bim").each(0,1,3) do |chr,snpname,position|
			keep << snpname unless @snp_particids[chr + ":" + position].nil?
		end

		File.open(@outputdir + "/keeplist.txt","w") do |keepfile|
			keep.each { |k| keepfile.puts k }
		end

		puts @snp_particids.keys.inspect
	end

	def recode_plink
		#Kernel.system("plink --noweb --bfile #{@plink_stem} --extract #{@outputdir}/keeplist.txt --recode12 --tab --out #{@outputdir}/barcode-snps")
	end

	def read_ped_genotypes 
		@snps = []
		DF.new(@outputdir + "/barcode-snps.map").each do |chr,snp,cm,pos| 
			@snps << chr.to_s + ":" + pos.to_s 
		end

		@particids = []
		@genotypes = []
		File.open(@outputdir + "/barcode-snps.ped") do |ped|
			ped.each do |line|
				fields = line.chomp.split("\t")
				@particids << fields.first
				6.times { fields.shift } 
				throw "Mismatch in number of genotypes #{fields.size} != #{@snps.size}" unless fields.size == @snps.size 
				@genotypes << fields.map { |f| nf = f.split(" ").sort.map { |g| g.to_i }.inject { |r,e| r+e }  }
			end
		end

		@particids_index = {}
		0.upto(@particids.size-1) { |i| @particids_index[@particids[i]] = i } 

		@snps_index = {}
		0.upto(@snps.size-1) { |i| @snps_index[@snps[i]] = i } 
		puts "-"*60
		puts @snps.inspect
		puts @snps_index.inspect

		@common_snps = @snps & @snp_particids.keys

		puts "Common SNPs: " + @common_snps.inspect

		# Compute frequency for each genotype:
		@snp_frequencies = {} 
		@common_snps.each do |snp|
			@snp_frequencies[snp] = {}
			snp_genotypes = []
			0.upto(@genotypes.size-1) do |i| 
				snp_genotypes << @genotypes[i][@snps_index[snp]]
			end
			snp_genotypes.uniq.each do |genotype|
					@snp_frequencies[snp][genotype] = (snp_genotypes.select { |g| g==genotype }).size.to_f / snp_genotypes.size.to_f
			end
		end
		puts "Common SNPs frequencies: " + @snp_frequencies.inspect

	end

	def make_numeric_geno 
		@snp_particids.keys.each do |snp|
			counts = {}
			@snp_particids[snp].each do |hsh|
				counts[hsh[:a1]] ||= 1
				counts[hsh[:a1]] = counts[hsh[:a1]] + 1 
				counts[hsh[:a2]] ||= 1
				counts[hsh[:a2]] = counts[hsh[:a2]] + 1 
			end

			# Determine major and minor allele
			major = (lambda { |cnts| cnts.inject { |p,n| n.last > p.last ? n : p }.first }).call(counts)
			minor = (lambda { |cnts| cnts.inject { |p,n| n.last < p.last ? n : p }.first }).call(counts)

			# Genotypes are coded using numerical values: major:major=3, major:minor=2, minor:minor=1
			@snp_particids[snp].each do |hsh|
				if hsh[:a1].nil? or hsh[:a2].nil? then
					hsh[:numeric_geno] = 0
				else
					hsh[:numeric_geno] = ((hsh[:a1] == major) ? 2 : 1) + ((hsh[:a2] == major) ? 2 : 1)
				end
			end
		end
	end

	def match_snp_geno
		# Match each SNP to its equivalent SNP in other dataset
		# This is to check that snp-wise correspondance is ok
		# If this is ok, we can then use to detect outlying individuals
		particid_mismatch = {}
		snp_stats = {}
		@common_snps.each do |snp|
			matching = 0
			mismatch = 0
			unknown = 0
			@snp_particids[snp].each do |hsh|
				#puts snp.inspect
				i,j = @particids_index[hsh[:particid]], @snps_index[snp]
				next if i.nil?

				if @genotypes[i][j].nil? or @genotypes[i][j] == 0 or hsh[:numeric_geno].nil? or hsh[:numeric_geno] == 0 then
					unknown = unknown + 1
				elsif @genotypes[i][j] == hsh[:numeric_geno] then
					matching = matching + 1
				else
					mismatch = mismatch + 1
					particid_mismatch[hsh[:particid]] ||= []
					particid_mismatch[hsh[:particid]] << snp
				end 
			end
			snp_stats[snp] = { :match => matching, :mismatch => mismatch, :missing => unknown, :total => matching+mismatch+unknown }
			puts snp_stats[snp].inspect
			puts "snp=#{snp} matching=#{matching} mismatch=#{mismatch} unknown=#{unknown}"
		end


		File.open(@output_file_mismatch,"w") do |file|
			file.puts "ID\tNMISMATCH\tERATE\tP-VALUE\tE-VALUE"
			(particid_mismatch.to_a.map { |k,v| [v.size,v,k] }).sort.reverse.each do |numsnps,snps,particid| 
				pvalue = 1
				# Calculate a p-value for particid (probability of seing this particular combination of mismatches at random) 
				snps.each { |snp| pvalue = pvalue * (snp_stats[snp][:mismatch].to_f / snp_stats[snp][:total].to_f) }
				# In the same breath, calculate average number of samples for the compared snps 
				avg_samples = (snps.map { |snp| snp_stats[snp][:total] }).inject { |r,v| r+v } / snps.size.to_f
				# Calculate e-value:
				evalue = avg_samples * pvalue
				file.puts "#{particid}\t#{numsnps}\t#{numsnps.to_f/@common_snps.size}\t#{pvalue}\t#{evalue}"
			end
		end
		puts "#{particid_mismatch.size} samples with barcoding errors written to #{@output_file_mismatch}"


		# Construct alignment matrix between all mismatching samples
		swaps = []
		particid_mismatch.keys.each do |particid_a|
			particid_mismatch.keys.each do |particid_b|
				next if particid_a == particid_b

				@particid_snps[particid_a] = @particid_snps[particid_a].delete_if { |x| x.nil? }

				unknown = 0
				matching = []
				mismatch = []
				@particid_snps[particid_a].each do |snp|
					#puts snp.inspect
					next unless @common_snps.include?(snp[:snpkey])
					i,j = @particids_index[particid_b], @snps_index[snp[:snpkey]] 
					if @genotypes[i][j].nil? or @genotypes[i][j] == 0 or snp[:numeric_geno].nil? or snp[:numeric_geno] == 0 then
						unknown = unknown + 1
					elsif @genotypes[i][j] == snp[:numeric_geno] then
						#matching << snp[:snpkey]
						matching << snp
					else
						mismatch << snp 
					end 
				end
				total = matching.size + unknown + mismatch.size 
				# Calculate p-value if we have a perfect match
				p_value = 1
				matching.each do |snp|
					#p_value = p_value * (snp_stats[snp][:match].to_f / snp_stats[snp][:total].to_f) 
					snpkey = snp[:snpkey]
					numeric_geno = snp[:numeric_geno]
					puts "#{snpkey} -- #{numeric_geno}"
					p_value = p_value * @snp_frequencies[snp[:snpkey]][snp[:numeric_geno]]
				end

				next if p_value >= 1
			
				adjust = 0

				# Monte Carlo simulation to calculate how often
				# we get a "match" p-value at least as extreme with as many matching snps
				if mismatch.size > 0 then
					binomials = []
					mismatch.size.downto(1) do |x|
						binomials << Binomial.new(mismatch.size + matching.size,x)
					end
					binomial_probs = normalize_cnts(binomials.map { |x| x.to_i } )
					smaller_probs = []
					larger_probs = []
					trials = 100*(1.to_f/p_value).to_i
					puts "trials: #{trials}"
					puts binomial_probs.inspect
					1.upto(trials) do |i|
						binom = binomials.random_select(binomial_probs)

						# Randomly select binom.bottom barcoding snps and calculate match probability	
						# record numbers of times that it is better or worse than matching p-value respective

						sample_probs = @particid_snps[particid_a].map { |x| @common_snps.include?(x[:snpkey]) ? (1.to_f / @common_snps.size) : 0 } 
						sampled_snps = @particid_snps[particid_a].sample_n(binom.bottom,sample_probs) 
					
						prob = sampled_snps.map do |x| 
							next if x.nil?
							if x.nil? then
								throw "x is nil" 
							end
							throw x if x[:snpkey].nil?  throw x if @snp_frequencies[x[:snpkey]].nil?
							y = @snp_frequencies[x[:snpkey]]
							throw y.to_s if y[x[:numeric_geno]].nil?
							y[x[:numeric_geno]]
						end

						# FIXME
						next if  prob.include?(nil)
			
						prob = prob.inject { |a,b| a*b }

						if prob > p_value then
							larger_probs << prob
						else
							smaller_probs << prob
						end
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

#				mismatch_p_value = 1
#				mismatch.each do |snp|
#					snpkey = snp[:snpkey]
#					numeric_geno = snp[:numeric_geno]
#					puts "#{snpkey} -- #{numeric_geno}"
#					mismatch_p_value = mismatch_p_value * (1-@snp_frequencies[snp[:snpkey]][snp[:numeric_geno]])
#				end

				# Adjust p-value by the number of mismatches
				#p_value = p_value * total.choose(matching.size)

#				combined_p_value = p_value / mismatch_p_value 

				#next if p_value > 1

				swaps << { :original_label => particid_a, :new_label => particid_b, :pvalue => p_value + adjust, :adjust => adjust, :evalue => particid_mismatch.size*p_value, :match => matching.size, :mismatch => mismatch.size } 
			end
		end
		File.open(@output_file_swaps, "w") do |file|
			keys = swaps.first.keys
			file.puts (keys.map { |x| x.upcase }).join("\t")
			swaps = swaps.sort { |a,b| a[:pvalue] <=> b[:pvalue] } 
			swaps.each { |swp| file.puts (keys.map { |k| swp[k] }).join("\t") }
		end
#		puts snp_stats.inspect
#		puts swaps.inspect
	end
end

barcoder = BarCoder.new(ARGV[0], ARGV[1], "output/barcoding")
#barcoder = BarCoder.new("HolbeckLGCCalls2.txt","update8","output/")
barcoder.load
barcoder.save_barcoding_data
barcoder.write_snp_keep_list
barcoder.recode_plink
barcoder.read_ped_genotypes
barcoder.make_numeric_geno
barcoder.match_snp_geno

