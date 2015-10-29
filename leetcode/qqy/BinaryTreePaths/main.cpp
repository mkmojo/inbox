class Solution {
public:
    vector<string> binaryTreePaths(TreeNode* root) {
        vector<string> v, vl, vr;
        if(root == NULL) {
            // do nothing to v
        }else {
            vl = binaryTreePaths(root->left);
            vr = binaryTreePaths(root->right);
            vl.insert(vl.end(), vr.begin(), vr.end());
            if(vl.empty()) v.push_back(to_string(root->val));
            else{
                for(auto &&it : vl){
                    string cnt_val = to_string(root->val);
                    string new_path = cnt_val + "->" + it;
                    v.push_back(new_path);
                }
            }
        }
        return v;
    }
};
