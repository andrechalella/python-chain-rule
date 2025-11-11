# Gemini - 2025-11-11

import importlib.util
import os
import sys

# 1. Get arguments and paths
if len(sys.argv) != 4:
    print("Error: Python script requires 3 arguments.", file=sys.stderr)
    sys.exit(1)

module_path = sys.argv[1]
temp_func_file = sys.argv[2]
temp_prog_file = sys.argv[3]

# 2. Dynamic Import
spec = importlib.util.spec_from_file_location("test_module", module_path)
if spec is None:
    print(f"Error: Could not find module at {module_path}", file=sys.stderr)
    sys.exit(1)

try:
    test_module = importlib.util.module_from_spec(spec)
    # The current directory is added to sys.path to resolve internal imports
    sys.path.insert(0, os.path.dirname(module_path) or '.')
    spec.loader.exec_module(test_module)
except Exception as e:
    print(f"Error loading module {module_path}: {e}", file=sys.stderr)
    sys.exit(1)

# 3. Extract symbols (lf and kwargs)
try:
    # lf: The LsodeFormatter instance
    lf = test_module.lf
    # make_program_file_kwargs: The dictionary of arguments
    make_program_file_kwargs = test_module.make_program_file_kwargs
except AttributeError as e:
    print(f"Error: Missing required symbol (lf or make_program_file_kwargs) in {module_path}. {e}", file=sys.stderr)
    sys.exit(1)

# 4. Run and output lf.make_functions_file()
output_func = lf.make_functions_file()
with open(temp_func_file, 'w') as f:
    f.write(output_func)

# 5. Run and output lf.make_program_file(**kwargs)
output_prog = lf.make_program_file(**make_program_file_kwargs)
with open(temp_prog_file, 'w') as f:
    f.write(output_prog)

# Clean up sys.path
sys.path.pop(0)
