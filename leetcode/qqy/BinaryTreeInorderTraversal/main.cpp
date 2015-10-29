class Solution {
    public:
        vector<int> inorderTraversal(TreeNode* root) {
            vector<int> ans;
            if(!root) return ans;
            stack<TreeNode*> pnodes;
            while( !pnodes.empty() || root != NULL){
                if(root){
                    pnodes.push(root);
                    root = root->left;
                }else{
                    root = pnodes.top();
                    pnodes.pop();
                    ans.push_back(root->val);
                    root = root->right;
                }
            }
            return ans;
        }
};
