#include<iostream>
#include<string>
#include<vector>
#include<set>
#include<map>
#include<algorithm>
using namespace std;

class Solution {

    int isPalindrome(string s){
        int i=0;
        int j=s.size();
        while(i<j){
            if(s[i] != s[j])
                return 0;
            i++;
            j--;
        }
        return s.size();
    }
    public:
    string longestPalindrome(string s) {
        if(s.size() <= 1) return s;
        int dp[s.size()][s.size()];
        for(size_t i=0;i<s.size();i++)
            for(size_t j=0;j<s.size();j++)
                dp[i][j] = -1;

        string res = "";
        for(size_t i=1;i<s.size();i++){
            for(size_t j=0;i-j+1>res.size();j++){
                dp[j][i] = 0;
                if(s[j] == s[i]){
                    if(dp[j+1][i-1] == -1){

                        dp[j][i] = (dp[j+1][i-1] = 
                                isPalindrome(s.substr(j+1, i-j-1))) + 2;

                        if( i-j+1 > res.size())
                            res = s.substr(j, i-j+1);
                    }else if(dp[j+1][i-1] == 0){
                        dp[j][i] = 0;
                    }else if(dp[j+1][i-1] > 0){
                        dp[j][i] = dp[j+1][i-1] + 2;
                        if( i-j+1 > res.size())
                            res = s.substr(j, i-j+1);
                    } else{
                        cout << "erroe" <<endl;
                    }


                }
            }
        }
        return res;
    }
};

int main()
{
    Solution s;
    cout << s.longestPalindrome("xxabccbatt") <<endl;
    cout << s.longestPalindrome("bccbatt") <<endl;
    cout << s.longestPalindrome("aaaaaaa") <<endl;
    cout << s.longestPalindrome("") <<endl;
    cout << s.longestPalindrome("xxabcdefgh") <<endl;
    cout << s.longestPalindrome("abcdefghxx") <<endl;
    cout << s.longestPalindrome("aaaacdefgbbbbbkkle") <<endl;
    return 0;
}

