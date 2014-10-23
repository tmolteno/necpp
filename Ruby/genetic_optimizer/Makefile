RUBY=ruby
OPTS=-I.
REXEC=${RUBY} ${OPTS}

#SERVER=localhost:8787
SERVER=`hostname`.`dnsdomainname`:8787
#SERVER=pugg.electron.otago.ac.nz:8787

.PHONY: test100 test500 test5000 test10000

test100:
	rm -rf test100
	mkdir -p test100
	${REXEC} Gp.rb --gen 10 --pop 100 --file test100/Gp

test500:
	rm -rf test500
	mkdir -p test500
	${REXEC} Gp.rb --clone 0.2 --mutate 0.2 --gen 100 --pop 500 --file test500/Gp --server ${SERVER}

test5000:
	rm -rf test5000
	mkdir -p test5000
	${REXEC} Gp.rb --clone 0.2 --mutate 0.2 --gen 200 --pop 2500 --file test5000/Gp --server ${SERVER}


test10000:
	rm -rf test10000
	mkdir -p test10000
	${REXEC} Gp.rb --clone 0.1 --mutate 0.2 --gen 500 --pop 50000 --file test10000/Gp --server ${SERVER}


client:
	${REXEC} SimulationClient.rb --server ${SERVER}

clean:
	- killall -9 ${RUBY}
	
rdoc:
	rm -rf doc
	rdoc --all --diagram --title 'RANT: Ruby Antenna Optimiser' README.md *.rb

unit_test:
	${REXEC} -d TestSuite.rb
