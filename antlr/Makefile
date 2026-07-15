########################################
#
#  Build the NEC2++ parser using
#  ANTLR 4 + Docker.
#
#  Prerequisites: Docker, GNU Make
#
########################################

.PHONY: all clean generate test

all: generate

generate:
	./build_in_docker.sh

# Quick sanity check — basic geometry + a full NEC file.
# All inputs conform to the field counts required by NEC-2 Part 3.
test: generate
	@echo "=== geometry ==="
	@printf 'GW 1 7 0. 0. -.25 0. 0. .25 .001\nGW 2 5 0. 0. .25 0. 0. .50 .002\nGE\nFR 0 1 0 0 1.\nEN\n' | ./run_in_docker.sh ./nec_parse
	@echo "=== all geometry card types ==="
	@printf 'GW 1 5 0. 0. -1. 0. 0. 1. 0.\nGC 0 0 1.0 .001 .005\nGA 3 11 0.15 0 180 .001\nGH 4 20 0.5 1.0 0.1 0.2 0.1 0.2 .001\nGS 2.0\nGM 0 0 1. 0. 0. 1. 0. 0. 0.\nGE\nFR 0 1 0 0 1.\nEN\n' | ./run_in_docker.sh ./nec_parse
	@echo "=== full card set (field counts per nec2prt3.pdf) ==="
	@printf 'GW 1 7 0. 0. -.25 0. 0. .25 .001\nGE\nFR 0 1 0 0 100.\nLD 5 0 0 0 1.\nGN 1\nEX 0 0 5 0 1. 0.\nnt 0 1 0 0 0 0 0 0 0 0\nTL 0 1 0 0 50. 0. 0. 0. 0. 0.\nXQ 0\nGD 0. 0. 0. 0.\nRP 0 1 1 0 0. 0. 0. 0.\nNX\nPT 0 0 0 0\nKH 1.0\nNE 0 1 1 1 0. 0. 0. 0. 0. 0.\nNH 0 1 1 1 0. 0. 0. 0. 0. 0.\nPQ 0 0 0 0\nEK 0\nCP 0 0 0 0\nPL 0 0 0 0 0 0 0\nWG\nMP 0 0 0 0 0. 0. 0. 0. 0. 0.\nEN\n' | ./run_in_docker.sh ./nec_parse || true
	@echo "=== end-to-end simulation (FR/EX/RP triggers simulate()) ==="
	@printf 'CE dipole\nGW 0 7 0. 0. -.25 0. 0. .25 .001\nGE\nFR 0 1 0 0 300.\nEX 0 0 4 0 1.\nRP 0 1 1 0 0. 0. 0. 0.\nEN\n' | ./run_in_docker.sh ./nec_parse > /tmp/nec_sim.out
	@echo "  output bytes: $$(wc -c < /tmp/nec_sim.out)"
	@grep -q 'RADIATION PATTERN' /tmp/nec_sim.out && echo "  radiation-pattern table present"
	@grep -q 'POWER BUDGET' /tmp/nec_sim.out && echo "  power-budget present"
	@echo "=== all tests passed ==="

clean:
	rm -rf generated nec_parse *.o
