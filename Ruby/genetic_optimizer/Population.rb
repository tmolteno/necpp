require("Individual.rb")
require 'thread'
require 'drb/drb'


class EvalQueue
	include DRb::DRbUndumped

	def initialize
		@q = Queue.new
	end

	def push(o)
		@q.push(o)
	end

	def empty?()
		@q.empty?()
	end

	def length()
		@q.length()
	end

	def pop
		@q.pop
	end
end

class IndividualArray < Array
    include Comparable 
end

class Population
	attr_accessor :individuals, :eval_queue

	include DRb::DRbUndumped

	# Server is a boolean
	def initialize(n, eval_queue, server=false)
		@n = n
		@individuals = IndividualArray.new
		@eval_queue = eval_queue
		@server = server
	end

	def randomize
		@n.times { |i|
			id = Individual.new
			id.randomize
			@individuals.push(id)
		}
		self.fitness
	end

	def random_fit_individual
		ret = nil
		sum = 0
		fmin = @individuals.last.fitness
		fs = @fitness_sum - @individuals.length*fmin
		r = rand(fs.floor)
#		print "fit sum = #{fs} fmin = #{fmin} rand = #{r}\n"
		@individuals.each { |ind|
			sum += ind.fitness - fmin
			# print "  fit sum(#{sorted[i].fitness - fmin}) = #{sum}\n"
			if (sum >= r)
				ret = ind.duplicate
				break
			end
		}
		ret = @individuals.last.duplicate if nil
		return ret
	end

    # Do the next generation
    def evolve(clone_fraction, mutant_fraction)
		print "evolve\n"

		next_gen = IndividualArray.new
		
		# Clone best fraction
		clones = (@n*clone_fraction).round
		print "Cloning #{clones} Individuals\n"

		0.upto(clones-1) { |i|
#			print "Clone #{i}\n"
			next_gen[i] = @individuals[i].duplicate
		}

		mutants = clones + (@n*mutant_fraction).round
		print "Mutating #{mutants - clones} Individuals #{clones}..#{mutants-1}\n"
		clones.upto(mutants-1) { |i|
			i1 = random_fit_individual
			i1.mutate
			next_gen[i] = i1
		}
		
		index = mutants
		print "Breeding #{mutants} .. #{@n}\n"
		# Breed the rest (selecting lower cost functions preferably)
		while (index < @n)
			i1 = random_fit_individual
			i2 = random_fit_individual
# 			print "Breed #{i1} #{i2}\n"
			c1,c2 = Individual.breed(i1,i2)
			next_gen[index] = c1
			next_gen[index+1] = c2 unless index >= @n-1
			index += 2
		end
		@individuals = next_gen
		print "Breeding Complete\n"
		self.fitness
		return @individuals[0].duplicate # return the best individual
	end
	
	def to_s
		ret = ""
		@individuals.each { |i| ret = ret + "#{i.fitness} #{i}" }
		return ret
	end

	def Population.sort(a, b)
		fa = a.fitness
		fb = b.fitness
		fa <=> fb
	end

	def fitness
		if (@server)
			# Add all the individuals to our queue (as fitness jobs)
			@individuals.each { |i|
				@eval_queue.push( i )
			}
			# Launch a thread to evaluate the fitnesses from the queue
			consumer = Thread.new do
				sleep_time = 3.0
				jobs_left = @eval_queue.length()
				jobs_per_second = -1.0
				until @eval_queue.empty? do
					print "Sleeping #{sleep_time}s: job_rate=#{sf(jobs_per_second,0.5)}, #{jobs_left} jobs left\n"
					sleep(sleep_time)
					jobs_per_second = (jobs_left - @eval_queue.length())/sleep_time
					jobs_left = @eval_queue.length()
					sleep_time = sf((jobs_left / jobs_per_second),0.2) unless jobs_per_second == 0
					sleep_time = 0.2 if sleep_time < 0.2
					sleep_time = 5.0 if sleep_time > 5.0
				end
			end
			consumer.join
			# Wait until the queue is empty
			print "Fitness Queue Done\n"
		else
			@individuals.each{ |a| a.fitness }
		end
		print "Sorting\n"
		@individuals.sort!
		@individuals.reverse!

		# Calculate the fitness sum (they should all be evaluated now)
		sum = 0
		@individuals.each { |i|
#			print "Fitness (#{i.fitness})\n"
			sum += i.fitness
		}
		print "Fitness Sum #{sum}\n"
		@fitness_sum = sum
	end
end