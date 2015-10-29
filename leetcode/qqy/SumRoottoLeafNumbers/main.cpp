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
    int str2num(string str){
        int res = 0;
        for(int i=0;i<str.size();i++){
            int cnt_val = str[i] - '0';
            res = res * 10 + cnt_val;
        }
        return res;
    }
public:
    int sumNumbers(TreeNode* root) {
        if(root == NULL) return 0;
        if(!root->left && !root->right) return root->val;
        string cnt_node = to_string(root->val);
        string left, right;
        if(!root->left && root->right) {
            left = "";
            right = to_string(sumNumbers(root->right));
            return str2num(cnt_node + right);
        }else if(root->left && !root->right){
            left = to_string(sumNumbers(root->left));
            right = "";
            return str2num(cnt_node + left);
        }else{
            left = to_string(sumNumbers(root->left));
            right = to_string(sumNumbers(root->right));
            return str2num(cnt_node + left) + str2num(cnt_node + right);
        }
       
    }
};

// I misunderstood the problem
// I thought each level will need to do the sum
// but instead, just need to do it at the top level


// The following is my AC version
// The idea is to fist find all possible paths then add them all up
// The problem is that the space consumption is high, can probabaly 
// do it partically while traversing. 
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
    int str2num(string str){
        int res = 0;
        for(int i=0;i<str.size();i++){
            int cnt_val = str[i] - '0';
            res = res * 10 + cnt_val;
        }
        return res;
    }
    vector<string> all_path(TreeNode* node){
        vector<string> ans;
        if(node == NULL) return ans;
        vector<string> left = all_path(node->left);
        vector<string> right = all_path(node->right);
        left.insert(left.end(), right.begin(), right.end());
        //need to be the first item if both children returns empty string
        if(left.empty()) ans.push_back(to_string(node->val));
        for(auto &&it:left){
            ans.push_back( to_string(node->val) + it);
        }
        return ans;
    }
public:
    int sumNumbers(TreeNode* root) {
        vector<string> paths = all_path(root);
        int ans = 0;
        for(auto &&it:paths){
            ans += str2num(it);
        }
        return ans;
    }
};
