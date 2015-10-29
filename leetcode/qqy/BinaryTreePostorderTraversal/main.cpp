class Solution {
    public:
        vector<int> postorderTraversal(TreeNode* root) {
            stack<TreeNode*> pn;
            vector<int> ans;
            TreeNode *last_vst; //this is useful when finished the right child
            if(root == NULL) return ans;
            while( !pn.empty() || root != NULL){
                if(root){
                    pn.push(root);
                    root = root->left;
                }else{
                    TreeNode *peek = pn.top();
                    if(peek->right != NULL && last_vst!=peek->right){
                        root = peek->right;
                    }else{
                        ans.push_back(peek->val);
                        last_vst = pn.top();
                        pn.pop();
                    }
                }
            }
            return ans;
        }
};
