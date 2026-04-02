/**
 * @file compact_defines.h
 * @author Zhuoyi Zou (zhuoyiz@andrew.cmu.edu)
 * 
 * Easier to read fixed-width types
 */

#ifndef _COMPACT_DEFINES_H_
#define _COMPACT_DEFINES_H_

#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <cstring>
#include <cstdlib>

/** @brief Inline min/max macros */
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/** @brief Inline rounding macro */
#define ROUNDUP(n,m) ((((n) + (m) - 1) / (m)) * (m))
#define ROUNDDOWN(n,m) (((n) / (m)) * (m))

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef int8_t i8;
typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;

#endif /* _COMPACT_DEFINES_H_ */
