PROGS   = hello pi_block_linear pi_block_tree pi_gather pi_nonblock_linear pi_one_side pi_reduce
all: $(PROGS)

%: %.cc
        mpicxx $< -o $@

clean:
        rm -f $(PROGS)