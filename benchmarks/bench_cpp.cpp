#include <benchmark/benchmark.h>

#include <cstdio>

void test_benchmark(benchmark::State& state)
{
	int a, b, c;
	for (auto _ : state)
	{
		c = a + b;
		benchmark::DoNotOptimize(c);
	}
}

BENCHMARK(test_benchmark);
