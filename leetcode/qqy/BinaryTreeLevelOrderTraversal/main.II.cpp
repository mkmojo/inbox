/**
 * Definition for a binary tree node.
 * struct TreeNode {
 *     int val;
 *     TreeNode *left;
 *     TreeNode *right;
 *     TreeNode(int x) : val(x), left(NULL), right(NULL) {}
 * };
 */
class Solution {
public:
    vector<vector<int>> levelOrderBottom(TreeNode* root) {
        vector<vector<int>> ans;
        if(root == NULL) return ans;
        vector<vector<int>> left = levelOrderBottom(root->left);
        vector<vector<int>> right = levelOrderBottom(root->right);
        //merge two lists from the bottom
        if(left.size() >= right.size()){
            int delta = left.size() - right.size();
            for(int i=0;i<delta;i++){
                ans.push_back(left[i]);
            }
            for(int i=0; i<right.size() ;i++){
                vector<int> entry_l = left[delta + i];
                vector<int> entry_r = right[i];
                entry_l.insert(entry_l.end(), entry_r.begin(), entry_r.end());
                ans.push_back(entry_l);
            }

        }else{
            int delta = right.size() - left.size();
            for(int i=0; i<delta; i++){
                ans.push_back(right[i]);
            }
            for(int i=0; i<left.size();i++){
                vector<int> entry_l = left[i];
                vector<int> entry_r = right[delta + i];
                entry_l.insert(entry_l.end(), entry_r.begin(), entry_r.end());
                ans.push_back(entry_l);
            }

        }

        ans.push_back(vector<int>({root->val}));
        return ans;
    }
};

//This first version has many problems:
//1. There is two copies of the code to do the merge depending on which 
//vector is larger in size.
//2. there is a lot of copying going on as the return parameter.

class Solution {
    vector<vector<int>> levelOrderBottom1(TreeNode* root) {
        vector<vector<int>> ans;
        if (root == NULL) return ans;
        ans.push_back(vector<int>({root->val}));
        vector<vector<int>> left = levelOrderBottom1(root->left);
        vector<vector<int>> right = levelOrderBottom1(root->right);
        //merge
        for(int i=0; i<left.size() || i<right.size();i++){
            if(i<left.size() && i<right.size()){
                vector<int> vel = left[i]; // left entry
                vector<int> ver = right[i];
                vel.insert(vel.end(), ver.begin(), ver.end());
                ans.push_back(vel);
            }else if(i < left.size()){
                ans.push_back(left[i]);
            }else{
                ans.push_back(right[i]);
            }
        }
        return ans;
    }
    public:
    vector<vector<int>> levelOrderBottom(TreeNode* root){
        vector<vector<int>> ans = levelOrderBottom1(root);
        reverse(ans.begin(), ans.end());
        return ans;
    }
};
//The thing worth learning is the way to organize the merge
//Notice that the condition for continuing the calculation
//is A OR B. Then in the else if arm, one of them is A. 
//
//This means else will handle all those satifies B. 
//
//This is very clean. 
