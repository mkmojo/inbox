#include<iostream>
#include<string>
#include<vector>
#include<set>
#include<map>
#include<algorithm>

using namespace std;
class Solution {
    string longestCommonPrefixForTwo(const string &a, const string &b){
        if(a == "" || b == "") return "";
        int n = a.size() > b.size() ? b.size() : a.size();
        int ans = 0;
        for(int i=n;i>=1;i--){
            string cntStr = a.substr(0,i);
            if(cntStr == b.substr(0,i)){
                ans = i;
                break;
            }
        }
        return a.substr(0, ans);
    }
public:
    string longestCommonPrefix(vector<string>& strs) {
        if(strs.size() == 0) return "";
        string ans = strs[0];
        for(int i=1;i<strs.size();i++){
            ans = longestCommonPrefixForTwo(ans, strs[i]);
            cout << "DEBUG: " << ans <<endl;
        }
        return ans;
    }
};

int main()
{
    Solution s;
    vector<string> x {"abddasdf", "abdda", "abd"};
    cout << s.longestCommonPrefix(x) <<endl;
    return 0;
}

