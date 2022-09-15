# gen single docs

SINGLE_F := RiskBased_intro.md Naive_intro.md Greedy_intro.md \
Util_readMkT.md Util_NYSEgen.md \

SINGLE_DOCS_out := $(addprefix $(OUTDIR)/, $(SINGLE_F))

define SINGLE_DOCS_make
$(eval mname := $(lastword $(subst /, ,$(1))))
$(eval dep := $(INDIR)/$(mname))

$(1): $(dep)
	@cat $(dep) > $(1)
endef
