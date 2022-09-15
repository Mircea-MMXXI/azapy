# gen risk-based docs

RB_F := CVaR MAD SMCR LSD MV SD GINI BTAD BTSD \

# rules RB_BACK_temp

RB_back_temp_Target := $(foreach ff, $(RB_F), $(TEMPDIR)/$(ff)_RB_back_temp.md)

define RB_back_temp_make
$(eval mname := $(firstword $(subst _, , $(lastword $(subst /, ,$(1))))))
$(eval dep := $(INDIR)/$(mname)_th_doc_base.md)

$(1): $(dep)
	@echo. > $(1)
	@echo [TOP]($(mname)_th_doc_base) >> $(1)
endef

# rule inline llink

RB_here_temp_Target := $(foreach ff, $(RB_F), $(TEMPDIR)/$(ff)_RB_here_temp.md)

define RB_here_temp_make
$(eval nname := $(firstword $(subst _, , $(lastword $(subst /, ,$(1))))))
$(eval dep := $(INDIR)/$(nname)Analyzer_class.md)

$(1): $(dep)
	@echo. > $(1)
	@echo Their meanings are [here]($(nname)Analyzer_class) >> $(1)
endef

# rules RB_RISK

RB_RISK := Risk_getWeights.md Risk_getRisk.md Risk_getPositions.md \
Risk_ViewFrontiers.md Risk_set_mktdata.md Risk_set_rrate.md \
Risk_set_rtype.md Risk_set_random_seed.md \

RB_RISK_temp_Target := $(foreach ff, $(RB_F), $(addprefix $(TEMPDIR)/_temp_$(ff)_, $(RB_RISK)))

define RB_RISK_temp_make
$(eval llink := $(firstword $(subst ., ,$(lastword $(subst _temp_, ,$(1))))))
$(eval mname := $(firstword $(subst _, ,$(llink))))
$(eval dep := $(INDIR)/$(lastword $(subst _$(mname)_, ,$(1))))

$(1): $(dep)
	@echo. > $(1)
	@echo ($(llink))= >> $(1)
endef

# rules RB_PORT

RB_PORT := Port_port_view.md Port_port_view_all.md Port_port_drawdown.md \
Port_port_perf.md Port_port_annual_returns.md Port_port_monthly_returns.md \
Port_port_period_returns.md Port_get_nshares.md Port_get_weights.md\
Port_get_account.md Port_get_mktdata.md \

RB_PORT_temp_Target := $(foreach ff, $(RB_F), $(addprefix $(TEMPDIR)/_temp_$(ff)_, $(RB_PORT)))

# rules RB_DOCS_RISK

RB_DOCS_RISK := _th_doc_base.md Analyzer_class.md Analyzer_class_example.md

RB_DOCS_RISK_temp_Target := $(foreach ff, $(RB_F), $(addprefix $(TEMPDIR)/_temp_$(ff),$(RB_DOCS_RISK)))

define RB_DOCS_RISK_temp_make
$(eval llink := $(firstword $(subst ., ,$(lastword $(subst _temp_, ,$(1))))))
$(eval dep := $(INDIR)/$(llink).md)

$(1): $(dep)
	@echo. > $(1)
	@echo ($(llink))= >> $(1)
endef

# rules RB_DOCS_PORT

RB_DOCS_PORT := _class.md _set_model.md _class_example.md

RB_DOCS_PORT_temp_Target := $(foreach ff, $(RB_F), $(addprefix $(TEMPDIR)/_temp_$(ff)_Port,$(RB_DOCS_PORT)))

define RB_DOCS_PORT_temp_make
$(eval llink := $(firstword $(subst ., ,$(lastword $(subst _temp_, ,$(1))))))
$(eval nname := $(firstword $(subst _, ,$(llink))))
$(eval fil := $(lastword $(subst _Port, ,$(llink))))
$(eval dep := $(INDIR)/Port_$(nname)$(fil).md)

$(1): $(dep)
	@echo. > $(1)
	@echo ($(llink))= >> $(1)
endef

# main risk-base docs

RB_DOCS_out := $(foreach ff, $(RB_F), $(OUTDIR)/$(ff)_th_doc.md)

define RB_DOCS_make
$(eval mname := $(firstword $(subst _, ,$(lastword $(subst /, ,$(1))))))

$(eval doclist := \
$(TEMPDIR)/_temp_$(mname)_th_doc_base.md $(INDIR)/$(mname)_th_doc_base.md $(TEMPDIR)/$(mname)_RB_back_temp.md \
$(TEMPDIR)/_temp_$(mname)Analyzer_class.md $(INDIR)/$(mname)Analyzer_class.md $(TEMPDIR)/$(mname)_RB_back_temp.md \
$(INDIR)/title_methods.md \
$(foreach ff, $(RB_RISK), $(TEMPDIR)/_temp_$(mname)_$(ff) $(INDIR)/$(ff) \
$(if $(filter Risk_getRisk.md, $(ff)), $(TEMPDIR)/$(mname)_RB_here_temp.md)  \
$(if $(filter Risk_getWeights.md, $(ff)), $(TEMPDIR)/$(mname)_RB_here_temp.md)  \
$(TEMPDIR)/$(mname)_RB_back_temp.md) \
$(TEMPDIR)/_temp_$(mname)Analyzer_class_example.md $(INDIR)/$(mname)Analyzer_class_example.md \
$(INDIR)/eb_prog.md $(ANALYZEREX)/$(mname)Analyzer_examples.py $(INDIR)/eb_prog.md $(TEMPDIR)/$(mname)_RB_back_temp.md \
$(TEMPDIR)/_temp_$(mname)_Port_class.md $(INDIR)/Port_$(mname)_class.md \
$(INDIR)/Port_constructor.md $(TEMPDIR)/$(mname)_RB_back_temp.md \
$(INDIR)/title_methods.md \
$(TEMPDIR)/_temp_$(mname)_Port_set_model.md $(INDIR)/Port_$(mname)_set_model.md $(TEMPDIR)/$(mname)_RB_back_temp.md \
$(foreach ff, $(RB_PORT), $(TEMPDIR)/_temp_$(mname)_$(ff) $(INDIR)/$(ff) $(TEMPDIR)/$(mname)_RB_back_temp.md) \
$(TEMPDIR)/_temp_$(mname)_Port_class_example.md $(INDIR)/Port_$(mname)_class_example.md \
$(INDIR)/eb_prog.md $(PORTFOLIOSEX)/Port_$(mname)_examples.py $(INDIR)/eb_prog.md $(TEMPDIR)/$(mname)_RB_back_temp.md)

$(1): $(sort $(doclist))
	@cat $(doclist) > $(1)
endef
