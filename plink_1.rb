# This class handles loading of plink files

class PlinkGenotypes
	attr_accessor
		:snps, 			# list of snps in plink-file (from .map file)
		:individuals, 		# list of individual ids (IID) in plink. We require these to be unique
		:genotypes,		# A two-dim array with @snps.size columns and @individuals.size rows 
		:snps_index,		# @snps_index maps from snp name to column id in @genotypes 
		:individuals_index 	# @individuals_index maps from IID to row index in @genotypes

	def init(plinkstem=nil,snps=nil) 
		@snps = []
		@individuals = [] 
		@genotypes = [] 
		@individuals_index = {}
		@snps_index = {}
		load12(plinkstem,snps) unless plinkstem.nil?
	end

	# Genotypes are expected to be in numeric format, as produced when using the --recode12
	# option of plink
	def load12(plinkstem,snps)  
		# Read plink .map to get a list of snps 
		DF.new(plinkstem + ".map").each do |chr,snp,cm,pos| 
			@snps << chr.to_s + ":" + pos.to_s 
		end

		# init @snps_index 
		0.upto(@snps.size-1) { |i| @snps_index[@snps[i]] = i } 

		# fill up genotypes and individuals from .ped file
		File.open(plinkstem + ".ped") do |ped|
			ped.each do |line|
				fields = line.chomp.split("\t")
				@individuals << fields.first
				6.times { fields.shift }  # Skip sample information fields
				throw "Mismatch in number of genotypes #{fields.size} != #{@snps.size}" unless fields.size == @snps.size 
				# Convert the two alleles to a number between 0-4, where 0=Missing both, 1=one allele missing, 2=AA, 3=AB, 4=BB
				# FIXME: one allele missing is crap (subtract one and abs val)
				@genotypes << fields.map { |f| nf = f.split(" ").sort.map { |g| g.to_i }.inject { |r,e| r+e }  }
			end
		end

		# init @individuals_index
		0.upto(@individuals_index.size-1) { |i| @individuals_index[@individuals[i]] = i } 


	end

	def common
		@common_snps = @snps & @snp_individuals.keys

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

end
