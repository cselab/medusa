#ifndef BOOLEAN_H
#define BOOLEAN_H

//ifdef __XLF
//define Boolean bool
//endif

// use this for a compiler which knows the bool type
//#define Boolean bool
// use this for a compiler which does not know it!
//ifndef __XLF
enum Boolean {FALSE=0, TRUE=1};
//endif
enum Sign {NEGATIVE = -1, ZERO = 0, POSITIVE = 1};

#endif


