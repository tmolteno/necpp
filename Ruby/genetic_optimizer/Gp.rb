require "Parameters"
require "Population"
require 'drb/drb'
require 'drb/acl'

# The methods that need to be overridden are
# Antenna.get_statistics(freq)
# * Individual.fitness_calc
# * Individual.fitness_from_stats(freq)
# TODO Fix up the fitness methonds above so that only one is overridden
class Gp
  include DRb::DRbUndumped

  #
  # Initialize the Genetic Programme
  #
  # * n - number of individuals in the starting population
  # * filename_base - the saved information from this run is stored in files that begin with this
  # * eval_queue - a queue that will be processed of individuals to be evaluated
  # * server - the name of the server that will manage the evaluations
  #
  def initialize(n, filename_base, eval_queue, server)
    @n = n

    @filename_base = filename_base
    @filename = "#{@filename_base}.dat"
    outfile = File.open(@filename, File::CREAT|File::TRUNC)
    outfile.close

    @population = Population.new(@n, eval_queue, server)
    @population.randomize
    @generation = 0
  end

  # Evolve the population for one generation, cloning a certain fraction
  # to appear in the next generation, and mutating a specified fraction, also
  # to appear in the next generation.
#
# This method works by telling the Population object to evolve
  def evolve(clone_fraction, mutant_fraction)
    @best = @population.evolve(clone_fraction, mutant_fraction)
    outfile = File.open(@filename, File::WRONLY|File::APPEND|File::CREAT) 
    outfile.print("#{@generation}, #{@best.fitness}\n")
    outfile.close
    geo = @best.get_antenna.get_nec_file(FREQUENCY)
    necFile = File.new("#{@filename_base}_#{@generation}.nec", "w")
    necFile.print(geo)
    necFile.close
    @generation += 1
  end

  def queue
    @population.eval_queue
  end
  
  def to_s
    @population.to_s
  end

  def best
    @best
  end
end

if __FILE__ == $0


  require 'optparse'
  options = {}
  options[:pop] = 100
  options[:gen] = 10
  options[:file] = "Gp_out"
  options[:clone] = 0.1
  options[:mutate] = 0.1
  begin
	  opts = OptionParser.new do |opts|
		  opts.banner = "Usage: Gp.rb [options]"

		  opts.on("--file S", String, "Use the specified file prefix") do |v|
			  options[:file] = v
		  end
		  opts.on("--pop N", Integer, "Population size") do |n|
			  options[:pop] = n
		  end
		  opts.on("--gen N", Integer, "Number of generations") do |n|
			  options[:gen] = n
		  end
		  opts.on("--clone N", Float, "Clone fraction [0..1]") do |n|
			  options[:clone] = n
		  end
		  opts.on("--mutate N", Float, "Mutation fraction [0..1]") do |n|
			  options[:mutate] = n
		  end
		  opts.on("--server S", String, "Run as a computation server (e.g. cyberiad.electron.otago.ac.nz:8787)") do |s|
			  options[:server] = s
		  end
	  end.parse!
  rescue Exception
	  puts "Error: Unknown command line option."
	  puts opts
	  exit
  end

  pop = options[:pop]
  gen = options[:gen]
  clone = options[:clone]
  mutate = options[:mutate]
  server = options[:server]

  print "RANT Genetic Programming Antenna Optimizer\n"
  print "Parameters pop=#{pop}, gen=#{gen}, clone=#{clone}, mutate=#{mutate}\n"
  
  eval_queue = EvalQueue.new


  if server != nil
	  URI="druby://#{server}"
	  acl = ACL.new(%w{deny all
	      allow localhost
	      allow 172.16.1/24})
	  #DRb.install_acl(acl)
	  print "Running as a compute server at #{URI}\n"
	  print "Run compute clients as follows:\n"
	  print "  ruby SimulationClient.rb --server #{server}\n"
	  DRb.start_service(URI, eval_queue)
  end
  g = Gp.new(pop, options[:file], eval_queue, (server != nil))
  
#	$SAFE = 1   # disable eval() and friends
  
  mainthread = Thread.new do
	  1.upto(gen) { |i|
		  print "generation ",  i, " \n"
		  g.evolve(clone, mutate)
		  print "\n\nBest Individual (#{g.best.fitness_display})\n#{g.best.get_antenna}\n"
	  }
  end
  # Wait for the drb server thread to finish before exiting.
  mainthread.join
  DRb.thread.exit unless server == nil
end
