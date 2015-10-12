#pragma once

#ifdef KORE_DEBUG
#define KORE_ASSERT(X, Y) if(X) {} else { throw (Y); }
#else 
#define KORE_ASSERT(X, Y) {}
#endif

namespace kore {
namespace debug {




}
} // namespace kore
