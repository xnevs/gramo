#ifndef EXPLORE_H_
#define EXPLORE_H_

template <
    typename State,
    typename Callback>
void explore(State & S, Callback callback = Callback()) {
  while(true) {
    bool S_full = S.full();
    if(S_full) {
      bool proceed = callback(S);
      if(!proceed) {
        break;
      }
    }
    if(!S_full && S.available()) {
      S.advance();
      bool success = S.assign();
      if(success) {
        S.push();
      } else {
        S.revert();
        S.next();
      }
    } else if(!S.empty()) {
      S.pop();
      S.revert();
      S.next();
    } else {
      break;
    }
  }
}

#endif  // EXPLORE_H_
