#include<iostream>
#include<string>
#include<vector>
#include<set>
#include<map>
#include<algorithm>

using namespace std;
class Solution {
    pair<int, int> grow(const string &s, int b){
        set<char> myset;
        int e = s.size()-1;
        for(int i=b;i<s.size();i++){
            myset.insert(s[i]);
            if(myset.size() != i-b+1){
                e = i - 1;
                break;
            }
        }
        //cout << "myset.size(): " <<myset.size() <<endl;
        //cout << "b e: " << b << e <<endl;
        return make_pair(b, e);
    }
public:    
    int len(const pair<int,int>  &p){
        return p.second - p.first + 1;
    }

    bool cmp(const pair<int,int> &p, const pair<int,int> &q){
        return len(p) < len(q);
    }

    int lengthOfLongestSubstring0(string s) {
        vector<pair<int, int> > vec;
        for(int i=0;i<s.size();i++){
            vec.push_back(grow(s, i));
        }
        sort(vec.begin(), vec.end(), 
                [](const pair<int,int> &p, const pair<int,int> &q){
                return (p.second - p.first) < (q.second - q.first); });

        for(int i=0;i<vec.size();i++){
            //cout << vec[i].first << " " << vec[i].second<<endl;
        }
        return len(vec[vec.size()-1]);
    }
////////////////////////////////////////////////////////////////////////////////

    int countDistinct(string s){
        set<char> myset;
        for(int i=0;i<s.size();i++){
            myset.insert(s[i]);
        }
        return myset.size();
    }

    int lengthOfLongestSubstring1(string s){
        int ans = 1;
        map<string, int> mymap;
        for(int i=0;i<s.size();i++){
            for(int j=s.size()-1;j>i;j--){
                string text = s.substr(i, j-i+1);
                if(mymap.find(text) == mymap.end()){
                    //have not seen this text before 
                    int nDistinct = countDistinct(text);
                    mymap[text] = nDistinct;
                }
                if(mymap[text] == text.size() && text.size() > ans){
                    ans = text.size();
                }
            }
        }
        return ans;
    }

    int lengthOfLongestSubstring2(string s){
        int ans = 1;
        map<string, int> mymap;
        for(int len=s.size();len>=1;len--){
            for(int i=0; i+len<=s.size();i++){
                string text = s.substr(i, len);
                //cout << text <<endl;
                if(mymap.find(text) == mymap.end()){
                    //have not seen this text before 
                    int nDistinct = countDistinct(text);
                    mymap[text] = nDistinct;
                }
                if(mymap[text] == text.size()){
                    return ans = len;
                }
            }
        }
        return ans;
    }

    int lengthOfLongestSubstring3(string s) {
        if(s.size() == 0) return 0;
        int dp = 1;
        map<string, int> mymap;
        for(int i=1;i<s.size();i++){
            for(int j = 0; i-j+1 > dp && j < i; j++){
                string text = s.substr(j, i-j+1);
                if(mymap.find(text) == mymap.end()){
                    int nDistinct = countDistinct(text);
                    mymap[text] = nDistinct;
                }
                if(mymap[text] == text.size() && text.size() > dp){
                    dp = text.size();
                }
            }
        }
        return dp;
    } 

    int lengthOfLongestSubstring(string s){
        int ans = 0;
        int b = 0;
        map<char, int> mymap;
        for(int i=0;i<s.size();i++){
            char c = s[i];
            int cntLen = 0;
            if(mymap.find(c) != mymap.end() && b <= mymap[c]){ 
                b = mymap[c] + 1;
            }

            cntLen = i - b + 1;
            if(cntLen > ans)
                ans = cntLen;

            mymap[c] = i; 
        }
        return ans;
    }
};

int main()
{
    Solution s;
    cout << s.lengthOfLongestSubstring("abcd") << endl;
    cout << s.lengthOfLongestSubstring("123456789") << endl;
    cout << s.lengthOfLongestSubstring("aaaa") << endl;
    cout << s.lengthOfLongestSubstring("lcumakzxkjhgyvbcbovuooifpybeeqdrsaetkfsvobdmwhqyvoujceot") << endl;
    cout << s.lengthOfLongestSubstring("vcfprwebjivnadykqplhzvmdjuttgsadwfsobyplgka") << endl;
    cout << s.lengthOfLongestSubstring("xbtowkgwssucmwqksalbxtplungmtklqlegkqmijhcnvafht") << endl;
    return 0;
}

