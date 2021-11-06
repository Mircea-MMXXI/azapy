
OUTDIR=../docs
INDIR=.

# gen risk-based docs

F_RB = CVaR MAD SMCR LSSD MV SD GINI SMGINI Omega

RISK_RB = Risk_getWeights.md Risk_getRisk.md Risk_getPositions.md \
Risk_ViewFrontiers.md Risk_set_mktdata.md Risk_set_rrate.md \
Risk_set_rtype.md Risk_set_random_seed.md

PORT_RB = Port_port_view.md Port_port_view_all.md Port_port_drawdown.md \
Port_port_perf.md Port_port_annual_returns.md Port_port_monthly_returns.md \
Port_port_period_returns.md Port_get_nshares.md \
Port_get_account.md Port_get_mktdata.md

define f_br_in_gen
	$(addprefix $(INDIR)/, $(1)_th_doc_base.md $(1)Analyzer_class.md $(RISK_RB) \
	Port_$(1)_class.md Port_$(1)_set_model.md $(PORT_RB))
endef

define f_out_gen
	$(OUTDIR)/$(1)_th_doc.md
endef

$(foreach ff, $(F_RB), $(eval $(ff)_in=$(call f_br_in_gen,$(ff))) )
$(foreach ff, $(F_RB), $(eval $(ff)_out=$(call f_out_gen,$(ff))) )

# gen Simple doc

F_SIMPLE = Simple

Port_Simple = Port_port_view.md Port_port_view_all.md Port_port_drawdown.md \
Port_port_perf.md Port_port_annual_returns.md Port_port_monthly_returns.md \
Port_get_mktdata.md

$(F_SIMPLE)_in = $(addprefix $(INDIR)/, Simple_th_doc_base.md Port_Simple_class.md \
Port_Simple_set_model.md $(Port_Simple))

$(F_SIMPLE)_out = $(call f_out_gen,Simple)

# gen naive docs

F_NAIVE = InvVol InvVar InvDD ConstW

PORT_NAIVE = Port_port_view.md Port_port_view_all.md Port_port_drawdown.md \
Port_port_perf.md Port_port_annual_returns.md  Port_port_monthly_returns.md \
Port_port_period_returns.md Port_get_nshares.md Port_get_account.md \
Port_get_mktdata.md

define f_naive_in_gen
	$(addprefix $(INDIR)/, $(1)_th_doc_base.md Port_$(1)_class.md \
	Port_$(1)_set_model.md $(PORT_NAIVE))
endef

$(foreach ff, $(F_NAIVE), $(eval $(ff)_in=$(call f_naive_in_gen,$(ff))) )
$(foreach ff, $(F_NAIVE), $(eval $(ff)_out=$(call f_out_gen,$(ff))) )

# gen singles

F_SINGLES = RiskBased_intro.md

$(F_SINGLES)_in = $(addprefix $(INDIR)/, $(F_SINGLES))
$(F_SINGLES)_out = $(addprefix $(OUTDIR)/, $(F_SINGLES))

# make targets

Target := $(F_RB) $(F_SIMPLE) $(F_NAIVE) $(F_SINGLES)

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
