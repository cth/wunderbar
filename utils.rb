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
		throw "Cannot samples #{n} from an array of size #{self.size}" if n > self.size
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


