pairs, size = [
  (2, 11),
  (3, 10),
  (6, 18),
  (7, 15),
  (8, 14),
  (12, 16),
  (13, 17),
], 18

def get_selection(seq_len, back, one_based_index):
    selection = set()
    def _solve(i, j):
        if (i, j) not in back:
            return
        t = back[i, j]
        if t == j:
            _solve(i, j - 1)
        else:
            offset = 1 if one_based_index else 0
            selection.add((t + offset, j + offset))
            _solve(i, t - 1)
            _solve(t + 1, j - 1)
    _solve(0, seq_len - 1)
    return selection

def get_db_structure(seq_len, pairs, one_based_index = True):
    pairs, selections = set(pairs), []
    offset = 1 if one_based_index else 0
    for itr in range(4):
      if len(pairs) == 0: # structure is complete
          break
      pairing = {}
      for pair in pairs:
          i, j = min(pair) - offset, max(pair) - offset
          pairing[i], pairing[j] = j, i
      
      dp = [[0 for _ in range(seq_len)] for _ in range(seq_len)]
      back = {}
      for span in range(2, seq_len + 1):
          for i in range(0, seq_len - span + 1):
              j = i + span - 1
              dp[i][j], back[i, j] = dp[i][j - 1], j # j unpaired
              t = pairing.get(j, -1)
              if t >= 0 and t >= i and t < j:
                  pair_reward = 1 + dp[i][t - 1] + dp[t + 1][j - 1] # j paired
                  if pair_reward > dp[i][j]:
                      back[i, j] = t
                      dp[i][j] = pair_reward
      
      selection = get_selection(seq_len, back, one_based_index)
      selections.append(selection)
      # print(f'page: {itr}, selection: {[(p, q) for p, q in selection]}')
      pairs = pairs - selection
    if len(pairs) > 0:
      print(f'[WARNING] Cannot incorporate all pairs in the structure, total pairs: {len(pairs)}, pairs left: {pairs}')
    
    db = ['.'] * seq_len
    for selection, symbol in zip(selections, [('(', ')'), ('[', ']'), ('{', '}'), ('<', '>')]):
      for i, j in selection:
          db[i - offset] = symbol[0]
          db[j - offset] = symbol[1]
    return ''.join(db)


if __name__ == '__main__':
  print(get_db_structure(size, pairs))
