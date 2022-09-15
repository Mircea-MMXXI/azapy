# gen GREEDY docs

GREEDY_F := Kelly \

# rules GREEDY_BACK_temp

GREEDY_back_temp_Target := $(foreach ff, $(GREEDY_F), $(TEMPDIR)/$(ff)_GREEDY_back_temp.md)

define GREEDY_back_temp_make
$(eval mname := $(firstword $(subst _, , $(lastword $(subst /, ,$(1))))))
$(eval dep := $(INDIR)/$(mname)_th_doc_base.md)

$(1): $(dep)
	@echo. > $(1)
	@echo [TOP]($(mname)_th_doc_base) >> $(1)
endef

# rules GREEDY_ENGINE

GREEDY_ENGINE :=

GREEDY_ENGINE_temp_Target := $(foreach ff, $(GREEDY_F), $(addprefix $(TEMPDIR)/_temp_$(ff)_, $(GREEDY_ENGINE)))

define GREEDY_ENGINE_temp_make
$(eval llink := $(firstword $(subst ., ,$(lastword $(subst _temp_, ,$(1))))))
$(eval mname := $(firstword $(subst _, ,$(llink))))
$(eval dep := $(INDIR)/$(lastword $(subst _$(mname)_, ,$(1))))

$(1): $(dep)
	@echo. > $(1)
	@echo ($(llink))= >> $(1)
endef

# rules GREEDY_PORT

GREEDY_PORT := Port_port_view.md Port_port_view_all.md Port_port_drawdown.md \
Port_port_perf.md Port_port_annual_returns.md Port_port_monthly_returns.md \
Port_port_period_returns.md Port_get_nshares.md Port_get_weights.md\
Port_get_account.md Port_get_mktdata.md \

GREEDY_PORT_temp_Target := $(foreach ff, $(GREEDY_F), $(addprefix $(TEMPDIR)/_temp_$(ff)_, $(GREEDY_PORT)))

# rules GREEDY_DOCS_RISK

GREEDY_DOCS_ENGINE := _th_doc_base.md Engine_class.md Engine_class_example.md

GREEDY_DOCS_ENGINE_temp_Target := $(foreach ff, $(GREEDY_F), $(addprefix $(TEMPDIR)/_temp_$(ff),$(GREEDY_DOCS_ENGINE)))

define GREEDY_DOCS_ENGINE_temp_make
$(eval llink := $(firstword $(subst ., ,$(lastword $(subst _temp_, ,$(1))))))
$(eval dep := $(INDIR)/$(llink).md)

$(1): $(dep)
	@echo. > $(1)
	@echo ($(llink))= >> $(1)
endef

# rules GREEDY_DOCS_PORT

GREEDY_DOCS_PORT := _class.md _set_model.md _class_example.md

GREEDY_DOCS_PORT_temp_Target := $(foreach ff, $(GREEDY_F), $(addprefix $(TEMPDIR)/_temp_$(ff)_Port,$(GREEDY_DOCS_PORT)))

define GREEDY_DOCS_PORT_temp_make
$(eval llink := $(firstword $(subst ., ,$(lastword $(subst _temp_, ,$(1))))))
$(eval nname := $(firstword $(subst _, ,$(llink))))
$(eval fil := $(lastword $(subst _Port, ,$(llink))))
$(eval dep := $(INDIR)/Port_$(nname)$(fil).md)

$(1): $(dep)
	@echo. > $(1)
	@echo ($(llink))= >> $(1)
endef

# main risk-base docs

GREEDY_DOCS_out := $(foreach ff, $(GREEDY_F), $(OUTDIR)/$(ff)_th_doc.md)

define GREEDY_DOCS_make
$(eval mname := $(firstword $(subst _, ,$(lastword $(subst /, ,$(1))))))

$(eval doclist := \
$(TEMPDIR)/_temp_$(mname)_th_doc_base.md $(INDIR)/$(mname)_th_doc_base.md $(TEMPDIR)/$(mname)_GREEDY_back_temp.md \
$(TEMPDIR)/_temp_$(mname)Engine_class.md $(INDIR)/$(mname)Engine_class.md \
$(TEMPDIR)/_temp_$(mname)Engine_class_example.md $(INDIR)/$(mname)Engine_class_example.md \
$(INDIR)/eb_prog.md $(ANALYZEREX)/$(mname)Engine_examples.py $(INDIR)/eb_prog.md $(TEMPDIR)/$(mname)_GREEDY_back_temp.md \
$(TEMPDIR)/_temp_$(mname)_Port_class.md $(INDIR)/Port_$(mname)_class.md \
$(TEMPDIR)/$(mname)_GREEDY_back_temp.md \
$(INDIR)/title_methods.md \
$(TEMPDIR)/_temp_$(mname)_Port_set_model.md $(INDIR)/Port_$(mname)_set_model.md $(TEMPDIR)/$(mname)_GREEDY_back_temp.md \
$(foreach ff, $(GREEDY_PORT), $(TEMPDIR)/_temp_$(mname)_$(ff) $(INDIR)/$(ff) $(TEMPDIR)/$(mname)_GREEDY_back_temp.md) \
$(TEMPDIR)/_temp_$(mname)_Port_class_example.md $(INDIR)/Port_$(mname)_class_example.md \
$(INDIR)/eb_prog.md $(PORTFOLIOSEX)/Port_$(mname)_examples.py $(INDIR)/eb_prog.md $(TEMPDIR)/$(mname)_GREEDY_back_temp.md)

$(1): $(sort $(doclist))
	@cat $(doclist) > $(1)
endef
