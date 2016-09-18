#include <array>
#include <algorithm>
#include <boost/sort/spreadsort/integer_sort.hpp>
#include "divsufsort_private.h"

static inline void sort_(int32_t* SA, int32_t* ISA, int32_t* PA, int64_t gidx, int64_t N) {
  //if (N < 2) __builtin_unreachable();
  auto* tmp = reinterpret_cast<std::tuple<int64_t, int32_t>*>(PA);

  for (int64_t i = 0; i < N; ++i) {
    int32_t a =  SA[i];
    int64_t b = ISA[a+1];
    int64_t c = ISA[a+2];
    std::get<0>(tmp[i]) = (b << 32) | c;
    std::get<1>(tmp[i]) = a;
  }
  boost::sort::spreadsort::integer_sort(tmp, tmp + N, [](auto a, int b){
    return std::get<0>(a) >> b;
  });

  // naming
  int64_t i = 0, j = gidx;
  auto t1 = tmp[0];
  auto t2 = tmp[0];
  for (; i < N-1; ++i) {
    ISA[std::get<1>(t1)] = j;
    // is the end of the group?
    t2 = tmp[i+1];
    if (std::get<0>(t1) != std::get<0>(t2)) { // upper bits unequal?
      j = gidx + 1 + i;
      std::get<1>(t1) = ~std::get<1>(t1);
    }
    SA[i] = std::get<1>(t1);
    t1 = t2;
  }
  SA[i] = ~std::get<1>(t1);
  ISA[~SA[i]] = j;
}

template<int M>
static inline void sort(int32_t* SA, int32_t* ISA, int64_t gidx, int64_t N) {
  //if (N < 2) __builtin_unreachable();
  //if (M < N) __builtin_unreachable();
  std::array<int64_t, M> tmp;

  if (M <= 16) {
    int64_t a =  SA[0];
    int64_t b = ISA[a+1];
    tmp[0] = a | (b << 32);
    for (int64_t i = 1, j; i < N; ++i) {
      a =  SA[i];
      b = ISA[a+1];
      auto v = a | (b << 32);  // direct insertion sort
      for (j = i; j >= 1 && v < tmp[j-1]; --j)
        tmp[j] = tmp[j-1];
      tmp[j] = v;
    }
  } else {
    for (int64_t i = 0; i < N; ++i) {
      int64_t a =  SA[i];
      int64_t b = ISA[a+1];
      tmp[i] = a | (b << 32);
    }
    std::sort(tmp.begin(), tmp.begin() + N);  // actually faster
  }

  // naming
  int64_t i = 0, j = gidx;
  auto t1 = tmp[0];
  auto t2 = tmp[0];
  for (; i < N-1; ++i) {
    ISA[static_cast<int32_t>(t1)] = j;
    // is the end of the group?
    t2 = tmp[i+1];
    if ((t1 ^ t2) >> 32) { // upper bits unequal?
      j = gidx + 1 + i;
      t1 = ~t1;
    }
    SA[i] = t1;
    t1 = t2;
  }
  SA[i] = ~t1;
  ISA[~SA[i]] = j;
}

void spsort(int32_t* ISA, int32_t* SA, int32_t* PA, int32_t n, int32_t m) {
  //printf("%i %i\n", n, m);
  for (int i = n - 1; 0 <= i; --i) {
    auto gidx   = ISA[i];     // get the group index (equals the start in SA)
    auto gfirst = SA + gidx;  // get the pointer to the start of the group

    if (*gfirst < 0) {        // check if already unique
      continue;               // sort next group
    }

    auto glast = gfirst;      // this could be skipped by rewriting the partition() function
    while (0 <= *++glast);    // search the end of the group
    *glast = ~*glast;         // flip to value to ease the rest of the code

    if (glast - gfirst + 1 <= 16) {
      sort<16>(gfirst, ISA, gidx, glast - gfirst + 1);
      continue;
    }

    if (glast - gfirst + 1 <= 1024) {
      sort<1024>(gfirst, ISA, gidx, glast - gfirst + 1);
      continue;
    }

    if (glast - gfirst + 1 <= m / 3) {
      sort_(gfirst, ISA, PA, gidx, glast - gfirst + 1);
      continue;
    }

    // this could be improved by switching to copying as soon as the chunks are small enough
    std::sort(gfirst, glast + 1, [&ISA](auto a, auto b) {
      if (ISA[a + 1] == ISA[b + 1]) return a < b; // this is needed to avoid naming problems
      return ISA[a + 1] < ISA[b + 1];
    });

    int64_t j = 0, k = gidx;
    auto t1 = gfirst[0];
    auto t2 = gfirst[0];
    for (; j < glast - gfirst; ++j) {
      t2 = gfirst[j+1];
      auto g1 = ISA[t1+1];
      auto g2 = ISA[t2+1];
      ISA[t1] = k;
      if (g1 != g2) {
        k = gidx + 1 + j;
        t1 = ~t1;
      }
      SA[j] = t1;
      t1 = t2;
    }
    SA[j] = ~t1;
    ISA[~SA[j]] = k;
  }
};
