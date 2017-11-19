#ifndef __HSA_H__
#define __HSA_H__

#include <cstring>
#include <iostream>
#include <cstdlib>
#include <queue>
#include <vector>
#include <math.h>
#include <fstream>
#include "SelectCase.h"
#include "NodeClass.h"
#include "BFTree.h"
#include "RandomTree.h"
#include "Setting.h"

double kCalculatePValue(double *numbers, int size, double testValue);
void kOutputAndGeneratePValueOne(std::string subRootS, std::string subRootT, char* pathCost, char* pathTreeS, char* pathTreeT, FILE* fp, int testNum);

namespace HSA {
    using TN::TreeNode;
    using TN::BFTree;
    using Setting::SCORE_MAX_G;
    using Setting::SCORE_MIN_G;
    using Setting::SCORE_MAX_L;
    using Setting::SCORE_MIN_L;
    
    const unsigned int LINE_MAX_C = 4096U;

    typedef Setting::CostType CostType;
    typedef Setting::ScoreType ScoreType;
    typedef SelectCase::BtType BtType;
    typedef NodeClass::Value NodeClassValue;
    typedef std::pair<size_t, size_t> Point;


	class HSA {
	//public:
		CostType costMatrix[Setting::CLASS_MAX][Setting::CLASS_MAX];
		BtType btMatrix[Setting::NODE_MAX][Setting::NODE_MAX];
		ScoreType scoreMatrix[Setting::NODE_MAX][Setting::NODE_MAX];
		NodeClass nodeClass;
		BFTree treeT;
		BFTree treeS;
		BFTree subTreeS, subTreeT;

		public:
			size_t TEST_CASE = 100U;

		inline bool BuildCost(char* path) {
			char buf[LINE_MAX_C];
			char aClass[LINE_MAX_C], bClass[LINE_MAX_C];
			ScoreType costScore;
			FILE* fp = fopen(path, "r");
			while (!feof(fp)) {
				if (fgets(buf, LINE_MAX_C, fp) &&
					sscanf(buf, "%s %s %d", aClass, bClass, &costScore) != 3) {
					fprintf(stderr, "Cost File Format Worry: %s\n", buf);
					continue;
				}
#ifdef DEBUG_Cost
				if (costScore < -128 || costScore > 127) {
					fprintf(stderr, "Cost File costScore Worry: %d\n", costScore);
					fclose(fp);
					return false;
				}
#endif // DEBUG_Cost
				NodeClassValue row = nodeClass.insert(std::string(aClass));
				NodeClassValue col = nodeClass.insert(std::string(bClass));
				costMatrix[col][row] = costMatrix[row][col] = costScore;
			}
			fclose(fp);

			/*std::cout << "\n costMatrix built is \n";
			for (int i = 0; i < 9; i++) {
				for (int j = 0; j < 9; j++) {
					std::cout << +costMatrix[i][j] << " ";
				}
				std::cout << "\n";
			}*/

			return true;
		}

		inline int BuildTree(char* path, BFTree& tree) {
#ifdef DEBUG_Tree
			int insertNum = 0;
#endif // DEBUG_Tree
			char buf[LINE_MAX_C];
			char lineage[LINE_MAX_C], name[LINE_MAX_C], _class[LINE_MAX_C];
			FILE* fp = fopen(path, "r");
			fgets(buf, LINE_MAX_C, fp);
			while (!feof(fp)) {
				fgets(buf, LINE_MAX_C, fp);
				if (sscanf(buf, "%s %s %s", lineage, name, _class) != 3) {
					fprintf(stderr, "Tree File Format Worry: %s\n", buf);
					fclose(fp);
					return -1;
				}
#ifdef DEBUG_Tree
				++insertNum;
#endif // DEBUG_Tree
				tree.insert(lineage, name, nodeClass.insert(std::string(_class)));
			}
			fclose(fp);
			tree.countAndSort();
#ifdef DEBUG_Tree
			if ((insertNum << 1) != (tree._size + 1)) {
				fprintf(stderr, "Tree File Worry: Not a BFTree.\n");
				return -1;
			}
			tree.print();
#endif // DEBUG_Tree
			return 0;
		}

		inline bool BuildTreeTS(char* pathS, char* pathT) {
			return BuildTree(pathS, treeS) == 0 && BuildTree(pathT, treeT) == 0;
		}

		inline void DP(void) {
			size_t m = (treeS._size + 1) >> 1, n = (treeT._size + 1) >> 1;
			TreeNode **pi = _row_, **pj = _col_;
			for (size_t i = 0; i < m; ++i, ++pi, pj = _col_)
				for (size_t j = 0; j < n; ++j, ++pj) {
					btMatrix[i][j] = SelectCase::NONE;
					scoreMatrix[i][j] = costMatrix[(*pi)->nodeClass][(*pj)->nodeClass];
				}
			pi = _row_; pj = _col_ + n;
			for (size_t i = 0; i < m; ++i, ++pi, pj = _col_ + n)
				for (size_t j = n; j < treeT._size; ++j, ++pj) {
					ScoreType lScore = _LSScore_;
					ScoreType rScore = _RSScore_;
					if (lScore > rScore) {
						btMatrix[i][j] = SelectCase::LS;
						scoreMatrix[i][j] = lScore;
					}
					else {
						btMatrix[i][j] = SelectCase::RS;
						scoreMatrix[i][j] = rScore;
					}
				}
			pi = _row_ + m; pj = _col_;
			for (size_t i = m; i < treeS._size; ++i, ++pi, pj = _col_)
				for (size_t j = 0; j < n; ++j, ++pj) {
					ScoreType lScore = _SLScore_;
					ScoreType rScore = _SRScore_;
					if (lScore > rScore) {
						btMatrix[i][j] = SelectCase::SL;
						scoreMatrix[i][j] = lScore;
					}
					else {
						btMatrix[i][j] = SelectCase::SR;
						scoreMatrix[i][j] = rScore;
					}
				}
			pi = _row_ + m; pj = _col_ + n;
			for (size_t i = m; i < treeS._size; ++i, ++pi, pj = _col_ + n)
				for (size_t j = n; j < treeT._size; ++j, ++pj) {
					BtType selectCase = SelectCase::NONE;
					ScoreType maxScore = SCORE_MIN_G;
					ScoreType score = _LSScore_;
					_Select_Max_(SelectCase::LS)
						score = _RSScore_;
					_Select_Max_(SelectCase::RS)
						score = _SLScore_;
					_Select_Max_(SelectCase::SL)
						score = _SRScore_;
					_Select_Max_(SelectCase::SR)
						score = _LLScore_;
					_Select_Max_(SelectCase::LL)
						score = _LRScore_;
					_Select_Max_(SelectCase::LR)
						_Select_Max_(SelectCase::SR)
						btMatrix[i][j] = selectCase;
					scoreMatrix[i][j] = maxScore;
				}
#ifdef DEBUG_Score
			static int db_score_time = 0;
			if (db_score_time++ > 0)
				return;
			for (size_t i = 0; i < treeS._size; ++i) {
				for (size_t j = 0; j < treeT._size; ++j) {
					printf("%d ", scoreMatrix[i][j]);
				}
				printf("\n");
			}
			printf("\n");
#endif // DEBUG_Score
		}

		inline ScoreType getMaxScore(void) {
			ScoreType maxScore = SCORE_MIN_G;
			for (size_t i = 0; i < treeS._size; ++i)
				for (size_t j = 0; j < treeT._size; ++j)
					if (scoreMatrix[i][j] > maxScore)
						maxScore = scoreMatrix[i][j];
			return maxScore;
		}

		inline void outputPValue(FILE* fp, ScoreType& maxScore) {
			fprintf(fp, "PValue:");
			double randomeTreeScores[Setting::MAX_TEST];//byK; Store the calculated scores
			ScoreType min = SCORE_MAX_G, max = SCORE_MIN_G, avg = 0;
			for (size_t i = 0; i < TEST_CASE; ++i) {
				TN::RandomTree tree(&treeT, nodeClass.size());
				tree.BuildTree();
				DP();
				ScoreType tmp = getMaxScore();
				if (tmp < min)
					min = tmp;
				if (tmp > max)
					max = tmp;
				avg += tmp;
				fprintf(fp, "%d ", tmp);
				randomeTreeScores[i] = tmp;//byK; Store the calculated scores 
			}

			fprintf(fp, "\nMin:%d Max:%d AVG:%d ", min, max, avg / TEST_CASE);
			//byK; calculate p-value
			if (avg / TEST_CASE > maxScore) {
				fprintf(fp, "SomethingWrong!\n");
			}
			else {
				double pvalue;
				pvalue = kCalculatePValue(randomeTreeScores, TEST_CASE, maxScore);
				fprintf(fp, "pvalue:%.5g \n", pvalue);
			}

		}


		inline void outputPValueLLL(FILE* fp, ScoreType& maxScore) {
			double randomeTreeScores[Setting::MAX_TEST];//byK; Store the calculated scores
			ScoreType min = SCORE_MAX_G, max = SCORE_MIN_G, avg = 0;
			fprintf(fp, "PValue:");
			for (size_t i = 0; i < TEST_CASE; ++i) {
				TN::RandomTree tree(&treeT, nodeClass.size());
				tree.BuildTree();
				DPL();
				ScoreType tmp = getMaxScore();
				if (tmp < min)
					min = tmp;
				if (tmp > max)
					max = tmp;
				avg += tmp;
				fprintf(fp, "%d ", tmp);
				randomeTreeScores[i] = tmp;//byK; Store the calculated scores 
			}

			fprintf(fp,"\nMin:%d Max:%d AVG:%d ", min, max, avg / TEST_CASE);
			//byK; calculate p-value
			double pvalue;
			pvalue = kCalculatePValue(randomeTreeScores, TEST_CASE, maxScore);
			fprintf(fp, "pvalue:%.5g \n", pvalue);
			

		}

		inline void getSubTree(std::string rootS, std::string rootT) {
			BFTree subTreeS, subTreeT;
			int sizeS = 0, sizeT = 0;
			TreeNode *subRankToNodeS[LINE_MAX_C], *subRankToNodeT[LINE_MAX_C];
			treeS.getRankToNodeWithID(subTreeS.rankToNode, rootS, &sizeS);
			treeT.getRankToNodeWithID(subTreeT.rankToNode, rootT, &sizeT);
			subTreeS._size = sizeS;
			subTreeT._size = sizeT;
			treeS = subTreeS;
			treeT = subTreeT;
			treeS.sort();
			treeT.sort();
			//std::ofstream myfile3;
			//myfile3.open("C:\\Users\\ATPs\\Documents\\GitHub\\VS\\HSA\\x64\\Release\\subTreeS.txt");
			//for (int i = 0; i < sizeS; i++) {
			//	std::cout << treeS.rankToNode[i]->id <<" ";
			//	//myfile3 << subRankToNodeS[i]->id <<"\n";
			//}
			//myfile3.close();
		}

	//	inline void outputPValueLT(std::string rootS, std::string rootT) {
	//		std::cout << rootS << " " << rootT << " test here!\n";
	//		double randomeTreeScores[TEST_CASE];//byK; Store the calculated scores
	//		ScoreType min = SCORE_MAX_G, max = SCORE_MIN_G, avg = 0;
	//		TreeNode *subRankToNodeS[LINE_MAX_C], *subRankToNodeT[LINE_MAX_C];
	//		int sizeS = 0, sizeT = 0;
	//		BFTree subTreeS, subTreeT;
	//		BFTree tmpTreeS, tmpTreeT;
	//		//BuildTree(pathS, tmpTreeS) == 0 && BuildTree(pathT, tmpTreeT) == 0
	//		treeS.getRankToNodeWithID(subRankToNodeS, rootS, &sizeS);
	//		std::cout << "\nsizeS is " << sizeS << "\n";
	//		treeS.getRankToNodeWithID(subRankToNodeT, rootT, &sizeT);
	//		std::cout << "\nsizeT is " << sizeT << "\n";
	//		
	//		//*subTreeS.rankToNode = *subRankToNodeS;
	//		subTreeS._size = sizeS;
	//		//*subTreeT.rankToNode = *subRankToNodeT;
	//		subTreeT._size = sizeT;
	//		std::cout << "size of subTreeS is " << subTreeS._size << "\n";
	//		std::cout << "size of subTreeT is " << subTreeT._size << "\n";


	//		for (int i = 0; i < sizeS; i++) {
	//			TreeNode *tempNode;
	//			tempNode = new TreeNode();
	//			//try to copy each of the TreeNode struct to avoid influence the original treeT and treeS
	//			tempNode->id = subRankToNodeS[i]->id;
	//			tempNode->name = subRankToNodeS[i]->name;
	//			tempNode->leaves = subRankToNodeS[i]->leaves;
	//			tempNode->nodeClass = subRankToNodeS[i]->nodeClass;
	//			tempNode->rank = subRankToNodeS[i]->rank;
	//			tempNode->left = subRankToNodeS[i]->left;
	//			tempNode->right = subRankToNodeS[i]->right;

	//			//tempNode = subRankToNodeS[i];
	//			subTreeS.rankToNode[i] = tempNode;
	//			//std::cout << subTreeS.rankToNode[i]->nodeClass <<" ";
	//		}
	//		subTreeS.sort();
	//		/*std::cout << "\n";
	//		for (int i = 0; i < sizeS; i++) {
	//			std::cout << subTreeS.rankToNode[i]->nodeClass << " ";
	//		}
	//		std::cout << "\n";*/

	//		for (int i = 0; i < sizeT; i++) {
	//			TreeNode *tempNode;
	//			tempNode = new TreeNode();
	//			//try to copy each of the TreeNode struct to avoid influence the original treeT and treeS
	//			tempNode->id = subRankToNodeT[i]->id;
	//			tempNode->name = subRankToNodeT[i]->name;
	//			tempNode->leaves = subRankToNodeT[i]->leaves;
	//			tempNode->nodeClass = subRankToNodeT[i]->nodeClass;
	//			tempNode->rank = subRankToNodeT[i]->rank;
	//			tempNode->left = subRankToNodeT[i]->left;
	//			tempNode->right = subRankToNodeT[i]->right;

	//			//tempNode = subRankToNodeS[i];
	//			subTreeT.rankToNode[i] = tempNode;
	//			//std::cout << subTreeT.rankToNode[i]->nodeClass << " ";
	//		}
	//		subTreeT.sort();
	//		/*std::cout << "\n";
	//		for (int i = 0; i < sizeT; i++) {
	//			std::cout << subTreeT.rankToNode[i]->nodeClass << " ";
	//		}*/

	//		for (int iii = 0; iii < 10; iii++) {
	//			/*std::ofstream myfile2;
	//			myfile2.open("C:\\Users\\ATPs\\Documents\\GitHub\\VS\\HSA\\x64\\Release\\treeT.txt");
	//			for (int i = 0; i < sizeT; i++) {
	//				myfile2 << subTreeT.rankToNode[i]->id << " ";
	//			}
	//			myfile2 << "\nafter randome\n";*/
	//			//randomize subTreeT
	//			BFTree *subTreeTloc;
	//			subTreeTloc = &subTreeT;
	//			TN::RandomTree subTree(subTreeTloc, nodeClass.size());
	//			subTree.BuildTree();
	//			/*for (int i = 0; i < sizeT; i++) {
	//				myfile2 << subTreeT.rankToNode[i]->id << " ";
	//			}
	//			myfile2.close();*/
	//			/*for (int i = 0; i < sizeT; i++) {
	//				if (std::cout << subTreeT.rankToNode[i]->name)
	//					std::cout << subTreeT.rankToNode[i]->nodeClass << " ";
	//			}*/

	//			//DPLrandTree(subTreeS, subTreeT, nodeClass, costMatrix);
	//			static BtType subBtMatrix[Setting::NODE_MAX][Setting::NODE_MAX];
	//			static ScoreType subScoreMatrix[Setting::NODE_MAX][Setting::NODE_MAX];


	//			//note: check the nodeClass size, composition, and costMatrix;
	//			std::cout << "\n nodeClass size is " << nodeClass.size() << "\n";
	//			////nodeClass seems good

	//			//check costMatrix;
	//			/*std::cout << "\n costMatrix is \n";
	//			for (int i = 0; i < 12; i++) {
	//			for (int j = 0; j < 12; j++) {
	//			std::cout << +costMatrix[i][j] << " ";
	//			}
	//			std::cout << "\n";
	//			}*/
	//			////costMatrix seems right. size of 9*9


	//			//dynamic programing
	//			size_t m = (subTreeS._size + 1) >> 1, n = (subTreeT._size + 1) >> 1;//checked, works fine.
	//			//std::cout << "size m " << m << " size n " << n << "\n";
	//			TreeNode **pi = _subrow_, **pj = _subcol_;
	//			///*std::cout << "\n nodeClass: ";
	//			//for (size_t i = 0; i < m; ++i, ++pi, pj = _subcol_) std::cout << (*pi)->nodeClass << " ";
	//			//std::cout << "\n nodeClass2: ";
	//			//for (size_t j = 0; j < n; ++j, ++pj) std::cout << (*pj)->nodeClass << " ";*/

	//			std::cout << "\n subScoreMatrix\n";
	//			for (int i = 0; i < m; ++i, ++pi, pj = _subcol_) {
	//				for (int j = 0; j < n; ++j, ++pj) {
	//					subBtMatrix[i][j] = SelectCase::NONE;
	//					int tempscore = costMatrix[(*pi)->nodeClass][(*pj)->nodeClass];
	//					//std::cout << tempscore << " ";
	//					//subScoreMatrix[i][j] = 1;//costMatrix[(*pi)->nodeClass][(*pj)->nodeClass];
	//					subScoreMatrix[i][j] = tempscore;//costMatrix[(*pi)->nodeClass][(*pj)->nodeClass];
	//				}
	//				//std::cout << "\n";
	//			}
	//			//delete subScoreMatrix;
	//			//delete subBtMatrix;

	//			pi = _subrow_; pj = _subcol_ + n;
	//			for (size_t i = 0; i < m; ++i, ++pi, pj = _subcol_ + n)
	//				for (size_t j = n; j < subTreeT._size; ++j, ++pj) {
	//					ScoreType lScore = _subLSScore_;
	//					ScoreType rScore = _subRSScore_;
	//					if (lScore > rScore) {
	//						if (lScore < 0) {
	//							subBtMatrix[i][j] = SelectCase::NONE;
	//							subScoreMatrix[i][j] = 0;
	//						}
	//						else {
	//							subBtMatrix[i][j] = SelectCase::LS;
	//							subScoreMatrix[i][j] = lScore;
	//						}
	//					}
	//					else {
	//						if (rScore < 0) {
	//							subBtMatrix[i][j] = SelectCase::NONE;
	//							subScoreMatrix[i][j] = 0;
	//						}
	//						else {
	//							subBtMatrix[i][j] = SelectCase::RS;
	//							subScoreMatrix[i][j] = rScore;
	//						}
	//					}
	//				}

	//			pi = _subrow_ + m; pj = _subcol_;
	//			for (size_t i = m; i < subTreeS._size; ++i, ++pi, pj = _subcol_)
	//				for (size_t j = 0; j < n; ++j, ++pj) {
	//					ScoreType lScore = _subSLScore_;
	//					ScoreType rScore = _subSRScore_;
	//					if (lScore > rScore) {
	//						if (lScore < 0) {
	//							subBtMatrix[i][j] = SelectCase::NONE;
	//							subScoreMatrix[i][j] = 0;
	//						}
	//						else {
	//							subBtMatrix[i][j] = SelectCase::SL;
	//							subScoreMatrix[i][j] = lScore;
	//						}
	//					}
	//					else {
	//						if (rScore < 0) {
	//							subBtMatrix[i][j] = SelectCase::NONE;
	//							subScoreMatrix[i][j] = 0;
	//						}
	//						else {
	//							subBtMatrix[i][j] = SelectCase::SR;
	//							subScoreMatrix[i][j] = rScore;
	//						}
	//					}
	//				}

	//			pi = _subrow_ + m; pj = _subcol_ + n;
	//			for (size_t i = m; i < subTreeS._size; ++i, ++pi, pj = _subcol_ + n)
	//				for (size_t j = n; j < subTreeT._size; ++j, ++pj) {
	//					BtType selectCase = SelectCase::NONE;
	//					ScoreType maxScore = SCORE_MIN_L;
	//					ScoreType score = _subLSScore_;
	//					_Select_Max_(SelectCase::LS)
	//						score = _subRSScore_;
	//					_Select_Max_(SelectCase::RS)
	//						score = _subSLScore_;
	//					_Select_Max_(SelectCase::SL)
	//						score = _subSRScore_;
	//					_Select_Max_(SelectCase::SR)
	//						score = _subLLScore_;
	//					_Select_Max_(SelectCase::LL)
	//						score = _subLRScore_;
	//					_Select_Max_(SelectCase::LR)
	//						subBtMatrix[i][j] = selectCase;
	//					subScoreMatrix[i][j] = maxScore;
	//				}
	//			std::ofstream myfile;
	//			myfile.open("C:\\Users\\ATPs\\Documents\\GitHub\\VS\\HSA\\x64\\Release\\costMatrix.tsv");
	//			std::cout << "\n subScoreMatrix: \n";
	//			ScoreType subMaxScoretmp, subMaxScore = SCORE_MIN_G;
	//			for (size_t i = 0; i < subTreeS._size; i++) {
	//				for (size_t j = 0; j < subTreeT._size; j++) {
	//					subMaxScoretmp = scoreMatrix[i][j];
	//					myfile << subMaxScoretmp << "\t";
	//					if (subMaxScoretmp > subMaxScore) {
	//						subMaxScore = scoreMatrix[i][j];
	//					}
	//				}
	//				myfile << "\n";
	//			}
	//			myfile.close();
	//			std::cout << "\n subMaxScore is " << subMaxScore << "\n";
	//		}
	//}




        inline void outputMatch(FILE* & fp, size_t & row, size_t & col, ScoreType& maxScore) {
            fprintf(fp, "Score:%d\nRootS:%s\nRootT:%s\n", maxScore,
                ((treeS.rankToNode[row])->id).c_str(),
                ((treeT.rankToNode[col])->id).c_str()
            );
            /*static_assert((sizeof(btMatrix[0]) == sizeof(btMatrix[0][0]) * Setting::NODE_MAX),
            "btMatrix Col NUM is not NODE_MAX");
            BtType* matrix = btMatrix[0];
            row = (row << Setting::NODE_MAX_OFFSET) | col;*/
            std::string pruneS, pruneT, matchS, matchT;
            std::queue<Point> que;
            do {
                switch (btMatrix[row][col]) {
                case SelectCase::NONE:
                    break;
                case SelectCase::LL: {
                    TreeNode* pi = _row_[row];
                    TreeNode* pj = _col_[col];
                    matchS += pi->left->id + " " + pi->right->id + " ";
                    matchT += pj->left->id + " " + pj->right->id + " ";
                    que.push(Point(pi->left->rank, pj->left->rank));
                    que.push(Point(pi->right->rank, pj->right->rank));
                    break;
                }
                case SelectCase::LR: {
                    TreeNode* pi = _row_[row];
                    TreeNode* pj = _col_[col];
                    matchS += pi->left->id + " " + pi->right->id + " ";
                    matchT += pj->right->id + " " + pj->left->id + " ";
                    que.push(Point(pi->left->rank, pj->right->rank));
                    que.push(Point(pi->right->rank, pj->left->rank));
                    break;
                }
                case SelectCase::LS: {
                    pruneT += _col_[col]->right->id + " ";
                    que.push(Point(_row_[row]->rank, _col_[col]->left->rank));
                    break;
                }
                case SelectCase::RS: {
                    pruneT += _col_[col]->left->id + " ";
                    que.push(Point(_row_[row]->rank, _col_[col]->right->rank));
                    break;
                }
                case SelectCase::SL: {
                    pruneS += _row_[row]->right->id + " ";
                    que.push(Point(_row_[row]->left->rank, _col_[col]->rank));
                    break;
                }
                case SelectCase::SR: {
                    pruneS += _row_[row]->left->id + " ";
                    que.push(Point(_row_[row]->right->rank, _col_[col]->rank));
                    break;
                }
                default:
                    fprintf(stderr, "btMatrix Worry In outputResult()\n");
                }
                if (que.empty()) {
                    break;
                } else {
                    std::pair<size_t, size_t> tmp = que.front();
                    row = tmp.first;
                    col = tmp.second;
                    que.pop();
                }
            } while (1);
            fprintf(fp, "PruneS:%s\nPruneT:%s\nMatchS:%s\nMatchT:%s\n",
                pruneS.c_str(), pruneT.c_str(), matchS.c_str(), matchT.c_str()
            );
        }

        inline void outputGResult(char* path) {
            DP();
            std::string outFile(path);
            FILE* fp = fopen((outFile + Setting::GlobalC).c_str(), "w");
            if (treeS._size == 0 || treeT._size == 0) {
                fprintf(stderr, "treeS._size == 0 || treeT._size == 0 In outputLResult()\n");
                return;
            }

            size_t row = treeS._size - 1, col = treeT._size - 1;
            ScoreType maxScore = scoreMatrix[row][col];
            
            outputMatch(fp, row, col, maxScore);
            outputPValue(fp,maxScore);
            fclose(fp);
        }

        inline std::vector<std::tuple<std::string, std::string>> outputLResult(char* path, size_t num) {
            std::string outFile(path);
			std::vector<std::tuple<std::string, std::string>> subTreeRootsTS;
            FILE* fp = fopen((outFile + Setting::LocalC).c_str(), "w");
            if (treeS._size == 0 || treeT._size == 0) {
                fprintf(stderr, "treeS._size == 0 || treeT._size == 0 In outputLResult()\n");
				return subTreeRootsTS;
            }

            DPL();
            size_t inum = 1U;
			ScoreType maxScoreMax;
			
			
        outputLResultForLoopEnd:
            while (inum <= num) {
#ifdef DEBUG_Bt
                for (size_t i = 0; i < treeS._size; ++i) {
                    for (size_t j = 0; j < treeT._size; ++j) {
                        printf("%d ", btMatrix[i][j]);
                    }
                    printf("\n");
                }
                printf("\n");
#endif // DEBUG_Bt
                size_t row, col;
                ScoreType maxScore = SCORE_MIN_L - 1;
                for (size_t i = 0; i < treeS._size; ++i)
                    for (size_t j = 0; j < treeT._size; ++j) {
                        if (_CanUseCase_(btMatrix[i][j]) && scoreMatrix[i][j] > maxScore) {
                            maxScore = scoreMatrix[i][j];
                            row = i;
                            col = j;
                        }
                    }
                if (maxScore < SCORE_MIN_L) {
                    fprintf(fp, "All node used.\n");
                    break;
                }

                std::string *rootS = &((treeS.rankToNode[row])->id);
                std::string *rootT = &((treeT.rankToNode[col])->id);
                std::string pruneS, pruneT, matchS, matchT;
                std::vector<Point> vec;
                vec.push_back(Point(row, col));
                size_t vecSize = 1U;
                do {

                    switch (_ReCase_(btMatrix[row][col])) {
                    case SelectCase::NONE:
                        break;
                    case SelectCase::LL: {
                        TreeNode* pi = _row_[row];
                        TreeNode* pj = _col_[col];
                        matchS += pi->left->id + " " + pi->right->id + " ";
                        matchT += pj->left->id + " " + pj->right->id + " ";
                        vec.push_back(Point(pi->left->rank, pj->left->rank));
                        vec.push_back(Point(pi->right->rank, pj->right->rank));
                        break;
                    }
                    case SelectCase::LR: {
                        TreeNode* pi = _row_[row];
                        TreeNode* pj = _col_[col];
                        matchS += pi->left->id + " " + pi->right->id + " ";
                        matchT += pj->right->id + " " + pj->left->id + " ";
                        vec.push_back(Point(pi->left->rank, pj->right->rank));
                        vec.push_back(Point(pi->right->rank, pj->left->rank));
                        break;
                    }
                    case SelectCase::LS: {
                        pruneT += _col_[col]->right->id + " ";
                        vec.push_back(Point(_row_[row]->rank, _col_[col]->left->rank));
                        break;
                    }
                    case SelectCase::RS: {
                        pruneT += _col_[col]->left->id + " ";
                        vec.push_back(Point(_row_[row]->rank, _col_[col]->right->rank));
                        break;
                    }
                    case SelectCase::SL: {
                        pruneS += _row_[row]->right->id + " ";
                        vec.push_back(Point(_row_[row]->left->rank, _col_[col]->rank));
                        break;
                    }
                    case SelectCase::SR: {
                        pruneS += _row_[row]->left->id + " ";
                        vec.push_back(Point(_row_[row]->right->rank, _col_[col]->rank));
                        break;
                    }
                    case SelectCase::USED: {
                        std::pair<size_t, size_t> tmp = vec.front();
                        btMatrix[tmp.first][tmp.second] += SelectCase::USED;
                        goto outputLResultForLoopEnd;
                    }
                    default:
                        fprintf(stderr, "btMatrix Worry In outputResult()\n");
                    }
                    if (vecSize == vec.size()) {
                        break;
                    } else {
                        std::pair<size_t, size_t> tmp = vec[vecSize++];
                        row = tmp.first;
                        col = tmp.second;
                    }
                } while (1);
                while (vecSize > 0) {
                    std::pair<size_t, size_t> tmp = vec[--vecSize];
                    btMatrix[tmp.first][tmp.second] = SelectCase::USED;
                }
                fprintf(fp,
                    "%u{\nScore:%d\nRootS:%s\nRootT:%s\nPruneS:%s\nPruneT:%s\nMatchS:%s\nMatchT:%s\n}\n",
                    inum, maxScore, rootS->c_str(), rootT->c_str(), 
                    pruneS.c_str(), pruneT.c_str(), matchS.c_str(), matchT.c_str()
                );


				maxScoreMax = maxScore;

				//outputPValueLT(rootS->c_str(), rootT->c_str());
				//std::cout << *rootS + " "+*rootT<<"\n";
				subTreeRootsTS.push_back(std::tuple<std::string, std::string>(*rootS, *rootT));


                ++inum;
            }
			
            
            fclose(fp);
			return subTreeRootsTS;
        }

        inline void DPL(void) {
                size_t m = (treeS._size + 1) >> 1, n = (treeT._size + 1) >> 1;
                TreeNode **pi = _row_, **pj = _col_;
                for (size_t i = 0; i < m; ++i, ++pi, pj = _col_)
                    for (size_t j = 0; j < n; ++j, ++pj) {
                        btMatrix[i][j] = SelectCase::NONE;
                        scoreMatrix[i][j] = costMatrix[(*pi)->nodeClass][(*pj)->nodeClass];
                    }
                pi = _row_; pj = _col_ + n;
                for (size_t i = 0; i < m; ++i, ++pi, pj = _col_ + n)
                    for (size_t j = n; j < treeT._size; ++j, ++pj) {
                        ScoreType lScore = _LSScore_;
                        ScoreType rScore = _RSScore_;
                        if (lScore > rScore) {
                            if (lScore < 0) {
                                btMatrix[i][j] = SelectCase::NONE;
                                scoreMatrix[i][j] = 0;
                            } else {
                                btMatrix[i][j] = SelectCase::LS;
                                scoreMatrix[i][j] = lScore;
                            }
                        } else {
                            if (rScore < 0) {
                                btMatrix[i][j] = SelectCase::NONE;
                                scoreMatrix[i][j] = 0;
                            } else {
                                btMatrix[i][j] = SelectCase::RS;
                                scoreMatrix[i][j] = rScore;
                            }
                        }
                    }
                pi = _row_ + m; pj = _col_;
                for (size_t i = m; i < treeS._size; ++i, ++pi, pj = _col_)
                    for (size_t j = 0; j < n; ++j, ++pj) {
                        ScoreType lScore = _SLScore_;
                        ScoreType rScore = _SRScore_;
                        if (lScore > rScore) {
                            if (lScore < 0) {
                                btMatrix[i][j] = SelectCase::NONE;
                                scoreMatrix[i][j] = 0;
                            } else {
                                btMatrix[i][j] = SelectCase::SL;
                                scoreMatrix[i][j] = lScore;
                            }

                        } else {
                            if (rScore < 0) {
                                btMatrix[i][j] = SelectCase::NONE;
                                scoreMatrix[i][j] = 0;
                            } else {
                                btMatrix[i][j] = SelectCase::SR;
                                scoreMatrix[i][j] = rScore;
                            }
                        }
                    }
                pi = _row_ + m; pj = _col_ + n;
                for (size_t i = m; i < treeS._size; ++i, ++pi, pj = _col_ + n)
                    for (size_t j = n; j < treeT._size; ++j, ++pj) {
                        BtType selectCase = SelectCase::NONE;
                        ScoreType maxScore = SCORE_MIN_L;
                        ScoreType score = _LSScore_;
                        _Select_Max_(SelectCase::LS)
                        score = _RSScore_;
                        _Select_Max_(SelectCase::RS)
                        score = _SLScore_;
                        _Select_Max_(SelectCase::SL)
                        score = _SRScore_;
                        _Select_Max_(SelectCase::SR)
                        score = _LLScore_;
                        _Select_Max_(SelectCase::LL)
                        score = _LRScore_;
                        _Select_Max_(SelectCase::LR)
                        btMatrix[i][j] = selectCase;
                        scoreMatrix[i][j] = maxScore;
                    }
#ifdef DEBUG_Score
                static int db_score_time = 0;
                if (db_score_time++ > 0)
                    return;
                for (size_t i = 0; i < treeS._size; ++i) {
                    for (size_t j = 0; j < treeT._size; ++j) {
                        printf("%d ", scoreMatrix[i][j]);
                    }
                    printf("\n");
                }
                printf("\n");
#endif // DEBUG_Score
        }

        HSA () {
            memset(costMatrix, 0, sizeof(costMatrix));//clean costMatrix
            for (size_t i = 0; i < Setting::CLASS_MAX; ++i)
                costMatrix[i][i] = Setting::MATCH_COST;
        }
        ~HSA () {
        }
    };
}

#endif // !__HSA_H__

double kCalculatePValue(double *numbers, int size, double testValue) {
	//with a list of numbers, size of numbers, and a testValue, return pvalue
	double pvalue, mean, sum = 0, sd = 0;
	for (int i = 0; i < size; i++) sum += numbers[i];
	mean = sum / size;
	for (int i = 0; i < size; i++) sd += pow(numbers[i] - mean, 2);
	sd = sqrt(sd / size);
	pvalue = erfc((testValue - mean) / sd);
	return pvalue;
}

void kOutputAndGeneratePValue(std::vector<std::tuple<std::string, std::string>> subTreeRootsST, char* pathCost, char* pathTreeS, char* pathTreeT, int testNum) {
	std::vector <std::tuple<std::string, std::string>> ::iterator it;
	int i = 1;
	std::string outFile(pathTreeS);
	FILE* fp = fopen((outFile + Setting::LocalC).c_str(), "a");
	for (it = subTreeRootsST.begin(); it<subTreeRootsST.end(); it++) {
		std::string subRootS, subRootT;
		subRootS = std::get<0>(*it);
		subRootT = std::get<1>(*it);
		fprintf(fp,"%d{\n",i);
		i++;
		kOutputAndGeneratePValueOne(subRootS, subRootT, pathCost, pathTreeS, pathTreeT, fp, testNum);
		fprintf(fp, "}\n");
	}
	fclose(fp);
}

void kOutputAndGeneratePValueOne(std::string subRootS, std::string subRootT, char* pathCost, char* pathTreeS, char* pathTreeT, FILE* fp, int testNum) {
	HSA::HSA *subHSA = new HSA::HSA();
	subHSA->TEST_CASE = testNum;
	subHSA->BuildCost(pathCost);
	subHSA->BuildTreeTS(pathTreeS, pathTreeT);
	subHSA->getSubTree(subRootS, subRootT);
	//std::string sa, sb;
	//std::cout << "subS: " << subRootS << " subT: " << subRootT << "\n";
	subHSA->DPL();
	int maxScore;
	maxScore = subHSA->getMaxScore();
	//std::cout << "\nmaxScore: " << maxScore << "\n";
	
	subHSA->outputPValueLLL(fp, maxScore);

	delete subHSA;
}
