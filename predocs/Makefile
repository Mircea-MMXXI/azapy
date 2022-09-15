OUTDIR=../docs
INDIR=.
TEMPDIR=../work/temp_docs
ANALYZEREX=../scripts/analyzers
PORTFOLIOSEX=../scripts/portfolios
UTILEX=../scripts/util

# rules

include Makefile_risk_base.mk Makefile_naive.mk Makefile_simple.mk
include Makefile_greedy.mk Makefile_single.mk Makefile_single_ex.mk

# build

all: $(RB_RISK_temp_Target) $(RB_back_temp_Target) $(RB_PORT_temp_Target) $(RB_DOCS_RISK_temp_Target) \
	$(RB_DOCS_PORT_temp_Target) $(RB_here_temp_Target) $(RB_DOCS_out) \
	$(NAIVE_back_temp_Target) $(NAIVE_PORT_temp_Target) $(NAIVE_DOCS_PORT_temp_Target) \
	$(NAIVE_DOCS_out) \
	$(SIMPLE_back_temp_Target) $(SIMPLE_PORT_temp_Target) $(SIMPLE_DOCS_PORT_temp_Target) \
	$(SIMPLE_DOCS_out) \
	$(GREEDY_ENGINE_temp_Target) $(GREEDY_back_temp_Target) $(GREEDY_PORT_temp_Target) \
	$(GREEDY_DOCS_ENGINE_temp_Target) $(GREEDY_DOCS_ENGINE_temp_Target) $(GREEDY_DOCS_out) \
	$(SINGLE_DOCS_out) \
	$(SINGLE_EX_temp_Target) $(SINGLE_EX_back_temp_Target) $(SINGLE_EX_DOCS_out)

$(foreach ff, $(RB_RISK_temp_Target), $(eval $(call RB_RISK_temp_make,$(ff))))
$(foreach ff, $(RB_back_temp_Target), $(eval $(call RB_back_temp_make,$(ff))))
$(foreach ff, $(RB_PORT_temp_Target), $(eval $(call RB_RISK_temp_make,$(ff))))
$(foreach ff, $(RB_DOCS_RISK_temp_Target), $(eval $(call RB_DOCS_RISK_temp_make,$(ff))))
$(foreach ff, $(RB_DOCS_PORT_temp_Target), $(eval $(call RB_DOCS_PORT_temp_make,$(ff))))
$(foreach ff, $(RB_here_temp_Target), $(eval $(call RB_here_temp_make,$(ff))))
$(foreach ff, $(RB_DOCS_out), $(eval $(call RB_DOCS_make,$(ff))))

$(foreach ff, $(NAIVE_back_temp_Target), $(eval $(call NAIVE_back_temp_make,$(ff))))
$(foreach ff, $(NAIVE_PORT_temp_Target), $(eval $(call NAIVE_PORT_temp_make,$(ff))))
$(foreach ff, $(NAIVE_DOCS_PORT_temp_Target), $(eval $(call NAIVE_DOCS_PORT_temp_make,$(ff))))
$(foreach ff, $(NAIVE_DOCS_out), $(eval $(call NAIVE_DOCS_make,$(ff))))

$(foreach ff, $(SIMPLE_back_temp_Target), $(eval $(call SIMPLE_back_temp_make,$(ff))))
$(foreach ff, $(SIMPLE_PORT_temp_Target), $(eval $(call SIMPLE_PORT_temp_make,$(ff))))
$(foreach ff, $(SIMPLE_DOCS_PORT_temp_Target), $(eval $(call SIMPLE_DOCS_PORT_temp_make,$(ff))))
$(foreach ff, $(SIMPLE_DOCS_out), $(eval $(call SIMPLE_DOCS_make,$(ff))))

#$(foreach ff, $(GREEDY_ENGINE_temp_Target), $(eval $(call GREEDY_ENGINE_temp_make,$(ff))))
$(foreach ff, $(GREEDY_back_temp_Target), $(eval $(call GREEDY_back_temp_make,$(ff))))
$(foreach ff, $(GREEDY_PORT_temp_Target), $(eval $(call GREEDY_ENGINE_temp_make,$(ff))))
$(foreach ff, $(GREEDY_DOCS_ENGINE_temp_Target), $(eval $(call GREEDY_DOCS_ENGINE_temp_make,$(ff))))
$(foreach ff, $(GREEDY_DOCS_PORT_temp_Target), $(eval $(call GREEDY_DOCS_PORT_temp_make,$(ff))))
$(foreach ff, $(GREEDY_DOCS_out), $(eval $(call GREEDY_DOCS_make,$(ff))))

$(foreach ff, $(SINGLE_DOCS_out), $(eval $(call SINGLE_DOCS_make,$(ff))))

$(foreach ff, $(SINGLE_EX_temp_Target), $(eval $(call SINGLE_EX_temp_make,$(ff))))
$(foreach ff, $(SINGLE_EX_back_temp_Target), $(eval $(call SINGLE_EX_back_temp_make,$(ff))))
$(foreach ff, $(SINGLE_EX_DOCS_out), $(eval $(call SINGLE_EX_DOCS_make,$(ff))))


.PHONY: clean cleandoc
clean:
	@rm -f $(RB_RISK_temp_Target)
	@rm -f $(RB_back_temp_Target)
	@rm -f $(RB_PORT_temp_Target)
	@rm -f $(RB_DOCS_RISK_temp_Target)
	@rm -f $(RB_DOCS_PORT_temp_Target)
	@rm -f $(RB_here_temp_Target)
	@rm -f $(RB_DOCS_out)


	@rm -f $(NAIVE_back_temp_Target)
	@rm -f $(NAIVE_PORT_temp_Target)
	@rm -f $(NAIVE_DOCS_PORT_temp_Target)
	@rm -f $(NAIVE_DOCS_out)

	@rm -f $(SIMPLE_back_temp_Target)
	@rm -f $(SIMPLE_PORT_temp_Target)
	@rm -f $(SIMPLE_DOCS_PORT_temp_Target)
	@rm -f $(SIMPLE_DOCS_out)

	@rm -f $(GREEDY_ENGINE_temp_Target)
	@rm -f $(GREEDY_back_temp_Target)
	@rm -f $(GREEDY_PORT_temp_Target)
	@rm -f $(GREEDY_DOCS_ENGINE_temp_Target)
	@rm -f $(GREEDY_DOCS_PORT_temp_Target)
	@rm -f $(GREEDY_DOCS_out)

	@rm -f $(SINGLE_DOCS_out)

	@rm -f $(SINGLE_EX_temp_Target)
	@rm -f $(SINGLE_EX_back_temp_Target)
	@rm -f $(SINGLE_EX_DOCS_out)


cleandoc:
	@rm -f $(RB_DOCS_out)
	@rm -f $(NAIVE_DOCS_out)
	@rm -f $(SIMPLE_DOCS_out)
	@rm -f $(GREEDY_DOCS_out)
	@rm -f $(SINGLE_DOCS_out)
	@rm -f $(SINGLE_EX_DOCS_out)
