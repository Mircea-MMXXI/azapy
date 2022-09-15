# gen SIMPLE docs

SIMPLE_F := Simple \

# rules SIMPLE_BACK_temp

SIMPLE_back_temp_Target := $(foreach ff, $(SIMPLE_F), $(TEMPDIR)/$(ff)_SIMPLE_back_temp.md)

define SIMPLE_back_temp_make
$(eval mname := $(firstword $(subst _, , $(lastword $(subst /, ,$(1))))))
$(eval dep := $(INDIR)/$(mname)_th_doc_base.md)

$(1): $(dep)
	@echo. > $(1)
	@echo [TOP]($(mname)_Port_th_doc_base) >> $(1)
endef

# rules SIMPLE_PORT

SIMPLE_PORT := Port_port_view.md Port_port_view_all.md Port_port_drawdown.md \
Port_port_perf.md Port_port_annual_returns.md Port_port_monthly_returns.md \
Port_get_mktdata.md \

SIMPLE_PORT_temp_Target := $(foreach ff, $(SIMPLE_F), $(addprefix $(TEMPDIR)/_temp_$(ff)_, $(SIMPLE_PORT)))

define SIMPLE_PORT_temp_make
$(eval llink := $(firstword $(subst ., ,$(lastword $(subst _temp_, ,$(1))))))
$(eval mname := $(firstword $(subst _, ,$(llink))))
$(eval dep := $(INDIR)/$(lastword $(subst _$(mname)_, ,$(1))))

$(1): $(dep)
	@echo. > $(1)
	@echo ($(llink))= >> $(1)
endef

# rules SIMPLE_DOCS_PORT

SIMPLE_DOCS_PORT := _th_doc_base.md _class.md _set_model.md _class_example_1.md _class_example_2.md

SIMPLE_DOCS_PORT_temp_Target := $(foreach ff, $(SIMPLE_F), $(addprefix $(TEMPDIR)/_temp_$(ff)_Port,$(SIMPLE_DOCS_PORT)))

define SIMPLE_DOCS_PORT_temp_make
$(eval llink := $(firstword $(subst ., ,$(lastword $(subst _temp_, ,$(1))))))
$(eval nname := $(firstword $(subst _, ,$(llink))))
$(eval fil := $(lastword $(subst _Port, ,$(llink))))
$(eval dep := $(if $(filter _th_doc_base, $(fil)), $(INDIR)/$(nname)$(fil).md, $(INDIR)/Port_$(nname)$(fil).md))

$(1): $(dep)
	@echo. > $(1)
	@echo ($(llink))= >> $(1)
endef

# main SIMPLE docs

SIMPLE_DOCS_out := $(foreach ff, $(SIMPLE_F), $(OUTDIR)/$(ff)_th_doc.md)

define SIMPLE_DOCS_make
$(eval mname := $(firstword $(subst _, ,$(lastword $(subst /, ,$(1))))))

$(eval doclist := \
$(TEMPDIR)/_temp_$(mname)_Port_th_doc_base.md $(INDIR)/$(mname)_th_doc_base.md $(TEMPDIR)/$(mname)_SIMPLE_back_temp.md \
$(TEMPDIR)/_temp_$(mname)_Port_class.md $(INDIR)/Port_$(mname)_class.md \
$(INDIR)/title_methods.md \
$(TEMPDIR)/_temp_$(mname)_Port_set_model.md $(INDIR)/Port_$(mname)_set_model.md $(TEMPDIR)/$(mname)_SIMPLE_back_temp.md \
$(foreach ff, $(SIMPLE_PORT), $(TEMPDIR)/_temp_$(mname)_$(ff) $(INDIR)/$(ff) $(TEMPDIR)/$(mname)_SIMPLE_back_temp.md) \
$(TEMPDIR)/_temp_$(mname)_Port_class_example_1.md $(INDIR)/Port_$(mname)_class_example_1.md \
$(INDIR)/eb_prog.md $(PORTFOLIOSEX)/Port_$(mname)_examples.py $(INDIR)/eb_prog.md $(TEMPDIR)/$(mname)_SIMPLE_back_temp.md \
$(TEMPDIR)/_temp_$(mname)_Port_class_example_2.md $(INDIR)/Port_$(mname)_class_example_2.md \
$(INDIR)/eb_prog.md $(PORTFOLIOSEX)/Port_$(mname)_examples2.py $(INDIR)/eb_prog.md $(TEMPDIR)/$(mname)_SIMPLE_back_temp.md)

$(1): $(sort $(doclist))
	@cat $(doclist) > $(1)
endef
