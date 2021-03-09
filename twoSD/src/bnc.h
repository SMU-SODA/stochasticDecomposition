/*
 * stoc.h
 *
 *  Created on: Sep 11, 2017
 *      Author: gjharsha
 */

#ifndef BNC_H_
#define BNC_H_

#include "twoSD.h"

#define maxdnodes   1000
#define useDNODE 
#define useINODE

#undef writeprob

int currKey;
int currDepth;
double GlobeUB;               // Global upper bound
oneProblem      *original;    // Info of the original problem 
struct BnCnodeType *bestNode; // best node that is found so far
struct BnCnodeType **nodearr; // array of deactivated leaf nodes
struct BnCnodeType **inodearr;// array of integer feasible leaf nodes
int dnodes;                   // number of deactivated nodes 
int inodes;                   // number of integer feasible nodes 
int maxcut;
double meanVal;               // Global lower bound 
#undef printBest
#undef printBranch
#undef depthtest
#undef writemaster

 /// BnC structures 
 /* A data structure which holds on the node informations in the branch and bound tree */
struct BnCnodeType {
	int		key;			             /* key id of the node. */
	int		parentkey;			         /* parent key id of the node. */
	int		depth;			             /* depth of the node on the tree. */
	int     varId;                       /* index of the variable that is disjuncted in this node */
	int     stInt;                       /* index of the first integer variable */
	int     edInt;                       /* index of the last integer variable */
	int     parentnumSamp;               /* sample size of the parent */
	int     numSamp;                     /* sample size of the current node */
	double  fracPi;                      /* (\tilde(Pi)_{n-1} + \Delta\Pi_{n})/\Pi_{n} */
	int     tightPi;                     /* number of pis used in tight cuts */
	int     Lambdasize;                  /* total lamda size */
	int     partightPi;                  /* number of pis used in tight cuts of th parent node */
	int     parLambdasize;               /* total lamda size the parent node */
	int     parparinit;                  /* initial lambda loop for the next iterations from parent of the parent */

	int     poolID;					 	 /* Index to cutsPool corresponding to the node. */
	int     parentPoolID;				 /* Index to cutsPool corresponding to the parent of the node. */

	iVector tcuts;						 /* collection of cuts indexes from cell->cuts are tight in this node */
	iVector partcuts;					 /* collection of parent cuts indexes from cell->cuts are tight in parent node */
	double  fracVal;                     /* disjuncted fractional value for the varId */
	double  parobjVal;                   /* objective value of the parent node */

	iVector  disjncs;                    /* list of disjunctive cuts on variables 0: not added 1: disjnct is added */
	dVector  * disjncsVal;               /* list of disjunctive cuts values on variables - it has upper and lower limits */

	iVector	IncumbiStar;				 /* indices of maximal pi for each distint observation for incumbent cuts */
	iVector	ParIncumbiStar;				 /* indices of maximal pi for each distint observation for incumbent cuts of the parent */
	int numVar;                          /* number of variables */
	int numRows;
	double LB;                           /* lower bounds and upper bounds at this node */
	double UB;                           /* lower bounds and upper bounds at this node */
	double objval;                       /* objective value of the node problem  */
	dVector vars;                        /* the values of variables */
//	dVector duals;                       /* the values of the duals of the node problem */
	bool   isInt;                        /* is the solution obtained from the node integer */
	bool   isActive;                     /* the node is active or not */
	bool   isSPopt;                      /* the estimation of the node completed */
	bool   isfathomed;                   /* the node is fathomed or not */
	bool   isleft;                       /* the node is in the left of its parent */
	bool   ishrstic;                     /* the node is a heuristic node */
	struct BnCnodeType *left, *right;    /* pointer to the left and right children. */
	struct BnCnodeType *nextnode;        /* next node in the queue. */
	struct BnCnodeType *prevnode;        /* next node in the queue. */
};


/* bnc.c */
int sumintVec(iVector a, int len);
bool isInVec(iVector vec, int len, int val);
struct BnCnodeType *newrootNode(int numVar, double LB, double UB, oneProblem * orig);
struct BnCnodeType *newNode(int key, struct BnCnodeType * parent, double fracVal, int varId, bool isleft);
int addBnCDisjnct(cellType *cell, dVector  *disjncsVal, int numCols, struct BnCnodeType * node);
int solveNode(stocType *stoc, probType **prob, cellType *cell, struct BnCnodeType *node, cString pname);
int setupNode(probType *prob, cellType *cell, struct BnCnodeType *node);
int cleanNode(probType *prob, cellType *cell, struct BnCnodeType *node);
void          freeNodes(struct BnCnodeType *root);
void 	     freeNode(struct BnCnodeType *node);
struct BnCnodeType *copyNode(struct BnCnodeType *node, double thresh);
struct BnCnodeType *nextNode(struct BnCnodeType *node);
struct BnCnodeType     *successorNode(struct BnCnodeType *node);
struct BnCnodeType *insertNode(struct BnCnodeType *node, struct BnCnodeType *activenode);
int branchbound(stocType *stoc, probType **prob, cellType *cell, double LB, double UB);
int branchVar(struct BnCnodeType *node, int strategy);
int branchNode(stocType *stoc, probType **prob, cellType *cell, struct BnCnodeType *node, struct BnCnodeType **activeNode);
int getfirstLeaf(int depth);
int getnodeIdx(int depth, int key, int isleft);
void fracLamda(cellType *cell, struct BnCnodeType *node);
void truncate(dVector var, dVector lb, dVector ub, int cnt);

#endif /* BNC_H_ */
