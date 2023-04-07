
#ifndef FASTCKY_UTILITY_V_H
#define FASTCKY_UTILITY_V_H

#include <string>
#include <vector>
#include <cmath>
#include <assert.h>

#define NTP(x,y) (x==1? (y==4?5:0) : (x==2? (y==3?1:0) : (x==3 ? (y==2?2:(y==4?3:0)) : (x==4 ? (y==3?4:(y==1?6:0)) : 0))))
#define PTLN(x) (x==1? 2:((x==2 || x==3)? 3:(x==5)? 1:4))
#define PTRN(x) (x==2? 2:((x==1 || x==4)? 3:(x==6)? 1:4))

#define NOTON 5 // NUM_OF_TYPE_OF_NUCS
#define NOTOND 25
#define NOTONT 125

#define EXPLICIT_MAX_LEN 4
#define SINGLE_MIN_LEN 0
#define SINGLE_MAX_LEN 20  // NOTE: *must* <= sizeof(char), otherwise modify State::TraceInfo accordingly

#define HAIRPIN_MAX_LEN 30
#define BULGE_MAX_LEN SINGLE_MAX_LEN
#define INTERNAL_MAX_LEN SINGLE_MAX_LEN
#define SYMMETRIC_MAX_LEN 15
#define ASYMMETRY_MAX_LEN 28
#define SPECIAL_HAIRPIN_SCORE_BASELINE -10000

extern bool _allowed_pairs[NOTON][NOTON];

#define MAXLOOP 30

#define GET_ACGU(x) ((x==1? 'A' : (x==2? 'C' : (x==3? 'G' : (x==4?'U': 'X')))))

#define GET_ACGU_NUC(x) ((x=='A'? 1 : (x=='C'? 2 : (x=='G'? 3 : (x=='U'?4: 0)))))

#define HAIRPINTYPE(x) ((x==5?0 : (x==6?1 : (x==8?2 : 3))))

extern int func1(std::string& a, int8_t b);

extern void func2(std::string& a, int b, std::vector<int>& c, std::vector<int>& d, std::vector<int>& e);

extern int func3(int a, int b, int c, int d, int e);

extern int func4(int a, int b, int c, int d, int e, int f, int g);

extern int func5(int a, int b, int c);

extern int func6(int a, int b, int c, int d, int e, int f, int g, int h);

extern int func7(int a, int b, int c, int d, int e, int h, int i);

extern int func8(int a, int b);

extern void func9(int a, int b);

extern int func10(int a, int b, int c);

extern int func11(int a, int b, int c, int d, int e, int f, int g, int h);

extern int func12(int a, int b, int c, int d, int e, int f, int g = -1);

extern int func13(int a, int b);

extern int func14(int a, int b, int c, int d, int e, int f, int g, int h, int i, int j, int k, int l);

extern int func15(int a, int b, int c, int d, int e, int f, int g);

#endif 

