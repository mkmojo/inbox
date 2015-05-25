#include<iostream>
#include<assert.h>
#include<string>
#include<algorithm>
#include<vector>
using namespace std;
class Solution {
    void swap(int&a, int&b){
        int t = a;
        a = b;
        b = t;
    }
public:
    vector<int> twoSum(vector<int>& nums, int target) {
        vector<int> ans;
        map<int, int> mymap;
        for(int i=0;i<nums.size();i++){
            int val = nums[i];
            int loc = i;
            mymap[val] = loc;
        }
        
        for(int i=0;i<nums.size();i++){
            int cnt = nums[i];
            int rest = target - cnt;
            if(mymap.find(rest) != mymap.end() && mymap[rest] != i){
                int small = i;
                int large = mymap[rest];
                if(small > large){
                    swap(small, large);
                }
                ans.push_back(small+1);
                ans.push_back(large+1);
                return ans;
            }
        }
    }
};

int main()
{
    Solution s;
    //Test for negative number
    //Test for small cases
    //Test for customized cases
    return 0; 
}
