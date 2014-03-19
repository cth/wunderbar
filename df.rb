# A simplistic sequential data file reader with some data.frame like capabilities,
# but without having to load the whole thing into memory
# Christian Theil Have

class DF
	def initialize(file,opts=nil) 
		@file = file
		@opts = opts.nil? ? {} : opts 
		@opts[:separator] ||= / +/
		@opts[:header] ||= false

		@field_names = []
		@field_name_indices = {}
	end

	def each(*fields)
		# automatically assume file with header if named fields are requested
		@opts[:header] ||= (fields.map { |f| f.class }).include?(String)

		File.open(@file) do |f|
			lineno = 0
			f.each do |line|
				line_fields = line.chomp.split(@separator).drop_while { |fn| fn =~ /^ *$/ }
				next if line_fields.empty?
				lineno = lineno + 1
				if lineno == 1 and @opts[:header] == true then
					@field_names = line_fields
					0.upto(@field_names.size-1) { |i| @field_name_indices[@field_names[i]] = i }
					next
				end
				#next if not @opts[:skip].nil? and lineno =< @opts[:skip] 

				extract_line_fields = []

				fields = 0.step(line_fields.size-1,1).to_a if fields.empty? 

				fields.each do |field|
					if field.class == Fixnum then
						extract_line_fields << line_fields[field]
					elsif field.class == String then
						throw "invalid field name '#{field}'" if @field_name_indices[field].nil?
						extract_line_fields << line_fields[@field_name_indices[field]]
					end
				end
				if extract_line_fields.length > 1
					yield extract_line_fields
				else
					yield extract_line_fields.first
				end
			end
		end
	end

	def select(fields) 

	end
end

# Example:
#DataFile.new("plink.lmiss").each(2,"SNP") do |chr,snp|
#	puts snp + " - " + chr
#end
