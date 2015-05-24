#include<iostream>
#include<string>
#include<algorithm>
#include<assert.h>
#include<vector>
#include<memory.h>
using namespace std;
class Solution {
    public:
        int numDecodings(string s) {
            if(s.size() <= 1) return s.size();
            vector<int> dp(s.size()+1, 0);
            dp[1] = 1;

            for(int i=2;i<=s.size();){
                int cntNum = 0;
                if(cntNum = stringToNumber(s.substr(i-1-1,2)) <= 26){
                    dp[i] = max(dp[i-2]+2, dp[i-1]+1);
                    i+=2;
                }

                dp[i] = dp[i-1]+1;
                i++;
            }


            return dp[s.size()-1];
        }

    private:
        int stringToNumber(string s){
            int sum = 0;
            for(int i=0;i<s.size();i++){
                int cntDig = charToNumber(s[i]);
                sum = 10*sum + cntDig;
            }
            return sum;
        }

        int charToNumber(char c){
            assert(c<='9' && c >= '0');
            return c - '0';
        }

        int printArray(int a[], int size){
            for(int i=0;i<size;i++){
                cout << a[i]<< " ";
            }
            cout << endl;
        }

};

int main()
{
    Solution s;
    cout<< s.numDecodings("12120")<< endl;
    return 0;
}

