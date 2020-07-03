#include <vector>
#include <iostream>
#include <string>


using namespace std;

vector<string> split_string(string);

// Complete the maxSubsetSum function below.

int maxSubsetSum(vector<int> arr) {
    if(arr.size() == 0) return 0;
    vector<int> dp(arr.size());
    dp[0] = arr[0];
    int ans = dp[0];

    if(ans < 2){
        return ans;
    }

    dp[1] = max(arr[1], arr[0]);
    ans = max(dp[1], ans);
    cout << "pni" << endl;
    for(int i = 2; i < arr.size(); i++){
        dp[i]=max(dp[i-2],max(dp[i-1],dp[i-2]+arr[i]));
        dp[i] = max(dp[i], arr[i]);
        cout << dp[i] << endl;
        ans = max(ans, dp[i]);
    }

    return ans;

}

int main()
{

  vector<int> arr;

  int n;

  cin >> n;

  int temp;
  for(int i = 0; i < n; i++){
    cin >> temp;

    arr.push_back(temp);
  }

cout << maxSubsetSum(arr) << endl;

    return 0;
}
