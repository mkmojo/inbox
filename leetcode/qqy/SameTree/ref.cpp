// This is a recursive solution
bool isSameTree(TreeNode *p, TreeNode *q){
    if(!p && !q) return true; 
    //Together with the last statement, this means one and only one tree is not empty. 
    if(!p || !q) return false; 
    return (p->val && q->val) &&
        isSameTree(p->left, q->left) &&
        isSameTree(p->right, q->right);
}

bool isSameTree(TreeNode *p, TreeNode *q){
    queue<TreeNode*>q1, q2;
    q1.push(p);
    q2.push(q);
    while( !q1.empty() && !q2.empty() ){
        TreeNode *cntp = q1.front();
        q1.pop();
        TreeNode *cntq = q2.front();
        q2.pop();
        if(!cntp && !cntq) continue;
        if(!cntp || !cntq) return false;

        if(cntp->val != cntq->val) return false;

        q1.push(cntp->left);
        q1.push(cntp->right);
        q2.push(cntq->left);
        q2.push(cntq->right);
    }
    return true;
}
