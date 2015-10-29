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
        vector<int> preorderTraversal(TreeNode* root) {
            stack<TreeNode*> pnodes;
            vector<int> ans;
            TreeNode* last_vst = NULL;
            if(root == NULL) return ans;
            while(!pnodes.empty() || root != NULL){
                if(root != NULL){
                    ans.push_back(root->val);
                    pnodes.push(root);
                    root = root->left;
                }else{
                    TreeNode *peek = pnodes.top();
                    if(peek->right != NULL && last_vst != peek->right)
                        root = peek->right;
                    else{
                        last_vst = pnodes.top();
                        pnodes.pop();
                    }
                }
            }
            return ans;
        }
};
