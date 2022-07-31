# gen single with examples

SINGLE_EX_F = Util_summary_MkTData.md Util_add_cash_security.md \
Util_update_all_MkTData.md \

SINGLE_EX_temp_Target := $(addprefix $(TEMPDIR)/_temp_,$(SINGLE_EX_F))

define SINGLE_EX_temp_make
$(eval llink := $(firstword $(subst ., ,$(lastword $(subst _temp_, ,$(1) )))))
$(eval dev := $(INDIR)/$(llink).md)

$(1): $(dev)
	@echo ($(llink)_TOP)= > $(1)
endef

SINGLE_EX_back_temp_Target := $(addprefix $(TEMPDIR)/_back_temp_,$(SINGLE_EX_F))

define SINGLE_EX_back_temp_make
$(eval llink := $(firstword $(subst ., ,$(lastword $(subst _temp_, ,$(1) )))))
$(eval dev := $(INDIR)/$(llink).md)

$(1): $(dev)
	@echo [TOP]($(llink)_TOP) > $(1)
endef

SINGLE_EX_DOCS_out := $(addprefix $(OUTDIR)/,$(SINGLE_EX_F))

define SINGLE_EX_DOCS_make
$(eval mname := $(firstword $(subst ., ,$(lastword $(subst Util_, ,$(1) )))))

$(eval doclist := $(TEMPDIR)/_temp_Util_$(mname).md $(INDIR)/Util_$(mname).md \
$(INDIR)/eb_prog.md $(UTILEX)/$(mname)_examples.py $(INDIR)/eb_prog.md \
$(TEMPDIR)/_back_temp_Util_$(mname).md)

$(1): $(sort $(doclist))
	@cat $(doclist) > $(1)
endef
