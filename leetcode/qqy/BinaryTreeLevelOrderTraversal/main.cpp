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
    vector<vector<int>> levelOrder(TreeNode* root) {
        vector<vector<int>> ans;
        if(root == NULL) return ans;
        queue<TreeNode*> q;
        q.push(root);
        TreeNode *lvl_last = root;
        vector<int> cnt_lvl;
        while(!q.empty()){
            TreeNode *cnt_node = q.front();
            q.pop();
            cnt_lvl.push_back(cnt_node->val);
            
            if(cnt_node->left) q.push(cnt_node->left);
            if(cnt_node->right) q.push(cnt_node->right);
            
            if(cnt_node == lvl_last){
                ans.push_back(cnt_lvl);
                cnt_lvl.clear();
                lvl_last = q.back();
            }
        }
        return ans;
    }
};
