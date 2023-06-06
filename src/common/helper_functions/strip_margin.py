def strip_margin(text: str) -> str:
  import re
  return re.sub("(^|\n)[ \t]*\|", "\\1", text)