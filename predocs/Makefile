
OUTDIR=../docs
INDIR=.

# gen risk-based docs

F_RB = CVaR MAD SMCR LSD MV SD GINI BTAD BTSD \

RISK_RB = Risk_getWeights.md Risk_getRisk.md Risk_getPositions.md \
Risk_ViewFrontiers.md Risk_set_mktdata.md Risk_set_rrate.md \
Risk_set_rtype.md Risk_set_random_seed.md \

PORT_RB = Port_port_view.md Port_port_view_all.md Port_port_drawdown.md \
Port_port_perf.md Port_port_annual_returns.md Port_port_monthly_returns.md \
Port_port_period_returns.md Port_get_nshares.md \
Port_get_account.md Port_get_mktdata.md \

ANALYZEREX=../scripts/analyzers/
PORTFOLIOSEX=../scripts/portfolios/

define f_br_in_gen
	$(addprefix $(INDIR)/, $(1)_th_doc_base.md $(1)Analyzer_class.md $(RISK_RB) \
	$(1)Analyzer_class_example.md $(ANALYZEREX)$(1)Analyzer_examples.py ex_1_end.md \
	Port_$(1)_class.md Port_constructor.md Port_$(1)_set_model.md $(PORT_RB) \
	Port_$(1)_class_example.md $(PORTFOLIOSEX)Port_$(1)_examples.py ex_2_end.md)
endef

define f_out_gen
	$(OUTDIR)/$(1)_th_doc.md
endef

$(foreach ff, $(F_RB), $(eval $(ff)_in=$(call f_br_in_gen,$(ff))) )
$(foreach ff, $(F_RB), $(eval $(ff)_out=$(call f_out_gen,$(ff))) )

# gen Simple doc

F_SIMPLE = Simple \

PORT_SIMPLE = Port_port_view.md Port_port_view_all.md Port_port_drawdown.md \
Port_port_perf.md Port_port_annual_returns.md Port_port_monthly_returns.md \
Port_get_mktdata.md \

define f_sp_in_gen
	$(addprefix $(INDIR)/, $(1)_th_doc_base.md \
	Port_$(1)_class.md Port_$(1)_set_model.md $(PORT_SIMPLE) \
	Port_$(1)_class_example_1.md \
	$(PORTFOLIOSEX)Port_$(1)_examples.py ex_2_end.md \
	Port_$(1)_class_example_2.md \
	$(PORTFOLIOSEX)Port_$(1)_examples2.py ex_2_end.md)
endef

$(foreach ff, $(F_SIMPLE), $(eval $(ff)_in=$(call f_sp_in_gen,$(ff))) )
$(foreach ff, $(F_SIMPLE), $(eval $(ff)_out=$(call f_out_gen,$(ff))) )

# gen naive docs

F_NAIVE = InvVol InvVar InvDD ConstW \

PORT_NAIVE = Port_port_view.md Port_port_view_all.md Port_port_drawdown.md \
Port_port_perf.md Port_port_annual_returns.md  Port_port_monthly_returns.md \
Port_port_period_returns.md Port_get_nshares.md Port_get_account.md \
Port_get_mktdata.md \

define f_naive_in_gen
	$(addprefix $(INDIR)/, $(1)_th_doc_base.md Port_$(1)_class.md \
	Port_constructor.md Port_$(1)_set_model.md \
	$(PORT_NAIVE) Port_$(1)_class_example.md \
	$(PORTFOLIOSEX)Port_$(1)_examples.py ex_2_end.md)
endef

$(foreach ff, $(F_NAIVE), $(eval $(ff)_in=$(call f_naive_in_gen,$(ff))) )
$(foreach ff, $(F_NAIVE), $(eval $(ff)_out=$(call f_out_gen,$(ff))) )

# gen singles

F_SINGLES = RiskBased_intro.md Naive_intro.md Greedy_intro.md \
Util_readMkT.md Util_NYSEgen.md \


$(foreach ff, $(F_SINGLES), $(eval $(ff)_in=$(INDIR)/$(ff)))
$(foreach ff, $(F_SINGLES), $(eval $(ff)_out=$(OUTDIR)/$(ff)))

# gen single with examples

UTILEX=../scripts/util/

F_SINGLES_EX = summary_MkTData add_cash_security \
update_all_MkTData \

define f_single_ex_in_gen
	$(addprefix $(INDIR)/, Util_$(1).md $(UTILEX)$(1)_examples.py ex_2_end.md)
endef

define f_single_ex_out_gen
	$(OUTDIR)/Util_$(1).md
endef

$(foreach ff, $(F_SINGLES_EX), $(eval $(ff)_in=$(call f_single_ex_in_gen,$(ff))) )
$(foreach ff, $(F_SINGLES_EX), $(eval $(ff)_out=$(call f_single_ex_out_gen,$(ff))) )

# gen greedy docs

F_GR = Kelly \

RISK_GR =

PORT_GR = Port_port_view.md Port_port_view_all.md Port_port_drawdown.md \
Port_port_perf.md Port_port_annual_returns.md Port_port_monthly_returns.md \
Port_port_period_returns.md Port_get_nshares.md \
Port_get_account.md Port_get_mktdata.md \

define f_gr_in_gen
	$(addprefix $(INDIR)/, $(1)_th_doc_base.md $(1)Engine_class.md $(RISK_GR) \
	$(1)Engine_class_example.md $(ANALYZEREX)$(1)Engine_examples.py ex_1_end.md \
	Port_$(1)_class.md Port_$(1)_set_model.md $(PORT_GR) \
	Port_$(1)_class_example.md $(PORTFOLIOSEX)Port_$(1)_examples.py ex_2_end.md)
endef

$(foreach ff, $(F_GR), $(eval $(ff)_in=$(call f_gr_in_gen,$(ff))) )
$(foreach ff, $(F_GR), $(eval $(ff)_out=$(call f_out_gen,$(ff))) )

# make targets

Target := $(F_RB), $(F_SIMPLE), $(F_NAIVE), $(F_SINGLES), $(F_SINGLES_EX), $(F_GR)

Target_out := $(foreach ff, $(addsuffix _out, $(Target)), $($(ff)))

define frule
$(1): $(2)
	@echo $(1)
	@cat $(2) > $(1)
endef

all:  $(Target_out)

$(foreach ff, $(Target), $(eval $(call frule, $($(ff)_out), $($(ff)_in))))

.PHONY: clean
clean:
	rm -f $(Target_out)
