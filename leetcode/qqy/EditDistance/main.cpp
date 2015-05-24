#include<iostream>
#include<algorithm>
#include<vector>
using namespace std;
class Solution {
    public:
        int minDistance(string word1, string word2) {
            if(word1.size() == 0 || word2.size()==0) 
                return max(word1.size(), word2.size());
            vector<vector<int> > dp(word1.size()+1, vector<int>(word2.size()+1, -1));
            //init matrix
            for(int i=0; i<=word1.size();i++){
                dp[i][0] = i;
            }

            for(int i=0; i<=word2.size();i++){
                dp[0][i] = i;
            }
            //perform dp
            for(int i=1; i<=word1.size();i++){
                for(int j=1;j<=word2.size();j++){
                    dp[i][j] = dp[i][j-1] + 1 < dp[i-1][j] + 1 ? 
                        min(dp[i-1][j-1] + (word1[i-1] == word2[j-1] ? 0 : 1), 
                                dp[i][j-1] + 1) :
                            min(dp[i-1][j-1] + (word1[i-1] == word2[j-1] ? 0 : 1), 
                                    dp[i-1][j] + 1);
                }
            }
            return dp[word1.size()][word2.size()];
        }
};

int main()
{
    Solution s;
    cout<< s.minDistance("pneumonoultramicroscopicsilicovolcanoconiosis","ultramicroscopic") << endl;
    cout<<endl;
    cout<< s.minDistance("abc","cbaasdfawefasdfasdf") << endl;
    cout<< s.minDistance("spartan","part") << endl;
}
