#
# The Action class represents actions for script specified on the command line
#
class Action
  def valid(action)
    self.actions.include? action.downcase.gsub("-", "_")
  end
  
  def actions
    (self.methods.map { |m| $1 if m =~ /^action_(.*)/ }).reject { |a| a.nil? }
  end
  
  def initialize(cmdline)
    @args = cmdline
    @actions, @options = [], {}
    
    begin
      raise "Insufficient arguments given!" if @args.length < 3
      @options[:datadir], @options[:outputfile] = @args.shift, @args.shift
      
      @args.each do |action|
        if action =~ /(.+)=(.+)/
          self.send("option_" + $1, $2)
        elsif valid(action) 
          @actions << "action_" + action
        else
          raise "!! Unknown action/option: #{action}"
        end end
      
      @reader = NCBIReader.new(@options)

    rescue => e
      puts "!" + e.to_s
      print_help
    end
  end
  
  def print_help
    File.open("README") do |f|
      f.each { |line| puts line } 
    end
    exit -1
  end


  def add_action(action)
    raise "#{action} is not a valid action" unless self.valid(action)
    @actions << "action_" + action.downcase.gsub("-", "_")
  end
  
  def run_actions
    @reader.process_data
    File.open(@options[:outputfile], "w") do |file|
      @actions.each do |a| 
        puts "Running #{a.gsub('_', ' ')}..."
        file << self.send(a).join("\n")
      end
    end
    puts "Wrote #{@options[:outputfile]}."
  end

  
