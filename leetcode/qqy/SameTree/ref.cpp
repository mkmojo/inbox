bool isSameTree(TreeNode *p, TreeNode *q){
    if(!p && !q) return true; 
    //Together with the last statement, this means one and only one tree is not empty. 
    if(!p || !q) return false; 
    return (p->val && q->val) &&
        isSameTree(p->left, q->left) &&
        isSameTree(p->right, q->right);
}
