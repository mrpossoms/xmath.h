.PHONY: format test

format:
	clang-format -i inc/*
	clang-tidy inc/* --checks=-*,clang-analyzer-core.*,clang-analyzer-cplusplus.*,clang-analyzer-nullability.*,clang-analyzer-valist.*

test:
	make -C tests test
