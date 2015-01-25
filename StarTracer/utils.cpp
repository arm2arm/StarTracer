
int rain () {
  CGalaxy *gal;
  int first[] = {-50,-10,15,20,25};
  int second[] = {-50,40,30,20,10};
  int nmax =sizeof(first)/sizeof(first[0]);
  vector<int> v(nmax);                           // 0  0  0  0  0  0  0  0  0  0
  vector<int>::iterator it;
  gal=new CGalaxy; 
  gal->SetIds(&first[0],nmax);
  int sgal=sizeof(gal);
  sgal=gal->id.size();
 
  sort (first,first+5);     //  5 10 15 20 25
  sort (second,second+5);   // 10 20 30 40 50

  it=set_intersection (first, first+5, second, second+5, v.begin());
                                               // 10 20 0  0  0  0  0  0  0  0

  cout << "intersection has " << int(it - v.begin()) << " elements.\n";
	
  it=v.begin();
  return 0;
}


int main(int, char **)
   {
   tree<string> tr;
   tree<string>::iterator top, one, two, loc, banana;

   tree<int> ht;
   tree<int>::iterator top_ht, one_ht, two_ht, loc_ht;

   top=tr.begin();
   one=tr.insert(top, "one");
   two=tr.append_child(one, "two");
   tr.append_child(two, "apple");
   banana=tr.append_child(two, "banana");
   tr.append_child(banana,"cherry");
   tr.append_child(two, "peach");
   tr.append_child(one,"three");
   loc= tr.append_child(one,"four");
   loc=tr.append_child(loc,"second four");
   loc=tr.append_child(loc,"third four");
   loc=tr.append_child(loc,"fourth four");
   loc=tr.append_child(loc,"fifth four");

   //////////////////////
   kptree::print_subtree_bracketed(tr, one);
   ////////////////////

   loc=find(tr.begin(), tr.end(), "one");
   if(loc!=tr.end()) {
      tree<string>::sibling_iterator sib=tr.begin(loc);
      while(sib!=tr.end(loc)) {
         cout << (*sib) << endl;
         ++sib;
         }
      cout << endl;
      tree<string>::iterator sib2=tr.begin(loc);
      tree<string>::iterator end2=tr.end(loc);
      while(sib2!=end2) {
         for(int i=0; i<tr.depth(sib2)-2; ++i) 
            cout << " ";
         cout << (*sib2) << endl;
         ++sib2;
         }
      }
   }

