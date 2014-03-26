# This class handles loading of plink files

class Genotypes
	attr_accessor :snps		# list of snps in plink-file (from .map file)
	attr_accessor :individuals	# list of individual ids (IID) in plink. We require these to be unique
	attr_accessor :genotypes	# A two-dim array with @snps.size columns and @individuals.size rows 
	attr_accessor :snps_index	# @snps_index maps from snp name to column id in @genotypes 
	attr_accessor :individuals_index # @individuals_index maps from IID to row index in @genotypes
	attr_accessor :snp_freq		# genotype frequencies by SNP

	def initialize(plinkstem=nil,snps=nil) 
		@snps = []
		@individuals = [] 
		@genotypes = [] 
		@individuals_index = {}
		@snps_index = {}
		load12(plinkstem,snps) unless plinkstem.nil?
		compute_frequencies
	end

	# Genotypes are expected to be in numeric, tab-separated format, as produced when using the --recode12 --tab
	# option of plink
	def load12(plinkstem,snps)  
		# Read plink .map to get a list of snps 
		DF.new(plinkstem + ".map").each { |chr,snp,cm,pos| @snps << chr.to_s + ":" + pos.to_s }

		# init @snps_index 
		0.upto(@snps.size-1) { |i| @snps_index[@snps[i]] = i } 

		rel_snp_ids = (snps.delete_if { |snp| @snps_index[snp].nil? }).map { |snp| @snps_index[snp] } unless snps.nil?
		
		# fill up genotypes and individuals from .ped file
		File.open(plinkstem + ".ped") do |ped|
			ped.each do |line|
				fields = line.chomp.split("\t")
				@individuals << fields[0] + " " + fields[1]
				6.times { fields.shift }  # Skip sample information fields
				throw "Mismatch in number of genotypes #{fields.size} != #{@snps.size}" unless fields.size == @snps.size
				# Select relevant fields
				if snps.nil? then
					relevant_fields = fields
				else
					relevant_fields = []
					0.upto(fields.size) { |i| relevant_fields << fields[i] if rel_snp_ids.include?(i) }
				end
	
				# Convert the two alleles to a number between 0-4, where 0=Missing, 1=AA, 2=AB, 3=BB
				@genotypes << relevant_fields.map { |f| nf = f.split(" ").sort.map { |g| g.to_i }.inject { |r,e| ((r+e)-1).abs }  }
			end
		end

		# init @individuals_index
		0.upto(@individuals.size-1) { |i| @individuals_index[@individuals[i]] = i } 

		# re-init @snps and @snps_index if necessary
		if not snps.nil? then
			new_snps = []
			new_snps_index = {}
			@snps.each do |snp| 
				new_snps << snp if snps.include?(snp)
				new_snps_index[snp] = new_snps.size-1
			end
			
			# overwrite old
			@snps_index = new_snps_index
			@snps = new_snps
		end
	end

	def compute_frequencies
		@snp_freq = {}
		@snps.each do |snp|
			@snp_freq[snp] = {}
			snp_id = @snps_index[snp]
			snp_genotypes = []
			0.upto(@genotypes.size-1) { |i| snp_genotypes << @genotypes[i][snp_id] }
			[0,1,2,3].each do |genotype|
					@snp_freq[snp][genotype] = (snp_genotypes.select { |g| g==genotype }).size.to_f / snp_genotypes.size.to_f
			end
		end
	end


	def each_by_snp(snp) 
		snp_id = @snps_index[snp]
		@individuals.each do |i|
			yield [ i, @genotypes[@individuals_index[i]][snp_id] ]
		end
	end

	def each_by_indv(indv)
		indv_id = @individuals_index[indv]
		@snps.each do |snp|
			snp_id = @snps_index[snp] 
			yield [ snp, @genotypes[indv_id][snp_id] ]	
		end
	end

	def include_individual?(indv)
		not @individuals_index[indv].nil?
	end

	def include_snp?(snp)
		not @snps_index[snp].nil?
	end

	def [](indv,snp) 
		#puts "indv: " + indv.to_s
		#puts "snp: " + snp.to_s
		genotypes[@individuals_index[indv]][@snps_index[snp]]
	end

	def to_s
		"#{individuals.size}x#{snps.size} genotype matrix"
		#"#{individuals.size}x#{snps.size} genotype matrix #{snps.inspect}"
	end
end
