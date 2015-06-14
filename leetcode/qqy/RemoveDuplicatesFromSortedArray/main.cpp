#include<iostream>
#include<string>
#include<algorithm>
#include<assert.h>
#include<vector>
#include<memory.h>
using namespace std;

class Solution {
    void printVec(vector<int>& v) const{
        for(int i=0;i<v.size();i++){
            cout << v[i] << " ";
        }
        cout <<endl;
    }

    void swap(vector<int>& v, int l, int r){
        int t = v[l];
        v[r] = t;
        v[l] = v[r];
    }

    public:
    int removeDuplicates1(vector<int>& nums) {
        if(nums.size() <= 1) return nums.size();
        int count = 1;
        int last = nums[0];
        int cnt = 0;
        vector<int>::iterator it;
        for(it = nums.begin(); it!=nums.end();it++){
            cnt = *it;
            if(cnt == last) nums.erase(it);
            last = cnt;
        }
        printVec(nums);
        return nums.size();
    }

    int removeDuplicates(vector<int>& nums) {
        if(nums.size() <= 1) return nums.size();

        int l = 0;
        int r = 0;
        while(r < nums.size()){
            if(nums[l] != nums[r]){
                swap(nums, l, r);
                l++;
                r++;
            }
            r++;
        }
        nums.resize(l+1);
        printVec(nums);
        return nums.size();
    }
};

int main()
{
    Solution s;
    vector<int> vec = {1, 1 , 2, 3 };
    vector<int> vec1 = {1, 1 , 2};
    vector<int> vec_empty;
    vector<int> vec_one = {1};
    vector<int> vec_all_one = {1,1,1,1,1,1,1,1,1,1};
    vector<int> vec_all_neg_one = {-1,  -1,  -1,  -1, -1};
    cout<< s.removeDuplicates(vec)<< endl;
    cout<< s.removeDuplicates(vec1)<< endl;
    cout<< s.removeDuplicates(vec_empty)<< endl;
    cout<< s.removeDuplicates(vec_one)<< endl;
    cout<< s.removeDuplicates(vec_all_one)<< endl;
    cout<< s.removeDuplicates(vec_all_neg_one)<< endl;
    return 0;
}

