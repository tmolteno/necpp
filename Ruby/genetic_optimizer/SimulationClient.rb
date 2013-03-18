require 'Population'
require 'Individual'
require 'drb/drb'
require 'Gp.rb'

class LineLogger
	def initialize
		@col = 0
		@tot = 0
	end

	def dot
		@tot += 1
		put('.')
	end

	def recount
		@tot = 0
		endl
	end

	def endl
		print ":#{@tot}\n"
		@col = 0
	end

	def put(text)
		@col += 1
		print "#{text}"
		endl if @col >= 75
		STDOUT.flush
	end
end

require 'optparse'

#
# An antenna simulation compute client, uses DRb to connect to a server
# and pops Individual antenna from the EvalQueue on the server, simulates
# them and returns an AntennaStatistics object
#
class SimulationClient

	def SimulationClient.run(server)
	    # The URI to connect to
#	    SERVER_URI="druby://#{server}"
	    server_uri="druby://#{server}"

	    # Start a local DRbServer to handle callbacks.
	    #
	    # Not necessary for this small example, but will be required
	    # as soon as we pass a non-marshallable object as an argument
	    # to a dRuby call.
	    DRb.start_service

	    log = LineLogger.new

	    while true do
	    begin
		    queue = DRbObject.new(nil,server_uri)
		    while true do
			    until queue.empty? do
				  begin
						individual = queue.pop
				    ant = individual.get_antenna
						s = Individual::get_stats(ant)
						individual.set_stats(s)
				    log.dot
					rescue Exception => e
						print "Exception: #{$!}\n"
						log.recount
						e.backtrace.each { |t|
						                   print "#{t}\n"
						                 }
						sleep(rand(10))
					end
			    end
		    end
	    rescue Exception => e
	    	print "Exception: #{$!}\n"
		    log.put("E.")
		    log.recount
	     	e.backtrace.each { |t|
	     		print "#{t}\n"
	     	}
		    sleep(rand(10))
	    end
	    end
	end
end

options = {}
begin
	opts = OptionParser.new do |opts|
		opts.banner = "Usage: SimulationClient.rb [options]"

		opts.on("--server S", String, "Connect to the computation server (e.g. cyberiad.electron.otago.ac.nz:8787)") do |s|
			options[:server] = s
		end
	end.parse!
rescue Exception
	    puts "Error: Unknown command line option."
	    puts opts
	    exit
end

server = options[:server]

unless server != nil
	puts "Error: A server URL must be provided."
	puts opts
	exit
end

def match_info(f, t)
	File.read(f).split("\n").each do |l|
	  if l.index(t) then
	    return l.split(":")[1].strip
		end
	end
	raise "ERROR: Cannot read `#{t}' from `#{f}'."
end

def cpu_cores
	begin
		cores = match_info('/proc/cpuinfo', 'cpu cores').to_i
	rescue
		cores = 1
	end
	cores
end

pids = Array.new
cpu_cores.times do
	job = fork do
		pids.push(SimulationClient.run(server))
	end
end

#pids.each { |pid| Process.wait(pid); }