import os
import sys

def find_F90_files():
    """
    The routine finds all the Fortran90 sources within
    the current working directory (cwd) and its subfolders.
    The current working directory is the location from which you run this scrip,
    typically the src directory of the project.
    The output is a dictionray (file_relativepath_dict) in which filenames are
    used as keys and the realtive path, w.r.t. the cwd, are values
    It uses the "os" package and its procdures to traversea all subdirectories.

    inputs: no input are necessary

    outputs:
        -file_relativepath_dict: the dictionray of filenames (keys) and
        relative paths (values)
    """

    file_relativepath_dict = dict()
    cwd = os.getcwd()
    for root, dirs, files in os.walk(cwd):
        for file in files:
            if file.endswith(".f90"):
                rel_path = os.path.relpath(root, cwd)
                file_relativepath_dict[file] = rel_path

    return file_relativepath_dict

def search_for_USE_and_MODULE_in_file(fname_with_path, exclude_use_list=list()):
    """
    This routine open the file named fname_with_path and searches inside it for
    used modules and decleared modules and submodules, and returns two sets:
    i) file_module_dependence, and ii) modules_declared_in_file.
    To generate such list this routine reads each line in the file
    and look for the USE, MODULE and SUBMODULE keywords. Furhtremore,
    some USE calls are associated to libraries, which usually are not in the
    main program directory, but are provided to the compiler via the -I and -L options.
    We want to discard libraries USE calls. All USE calls we want to must be
    provided via the exclude_use_list. For instance exclude_use_list=[omp_lib] will
    skip all USE omp_lib statements.


    inputs:
        - fname_with_path: the name of the file, including its path, that this routine parses

        - exclude_use_list (optional): the list of modules we want to skip while tracking USE calls.

    outputs:
        -file_module_dependence: the set of modules used by this file

        -modules_declared_in_file: the set of modules and submodules decleared by this file
    """
    # intiialize the set containing the list of submodule each files depends on
    file_module_dependence = set()
    modules_declared_in_file = set()
    module_submodule_dependence = set()
    submodules_declared_in_file = set()

    # when we look for a module declaration we usally look for
    # module xxx, however there a few keywords we want to avoid.
    # We include all such keywords in the following list
    module_exclude_keywords= ["procedure", "subroutine", "function"]

    # process the file line by line
    with open(fname_with_path,'r') as fin:
        for line in fin:
            # first clean the line by removing unwanted character, this might also be
            # done via regular expression but for now we are lazy
            line=line.replace('\t', ' ')
            line=line.replace(',', ' ')
            line=line.replace(':', ' ')
            line=line.replace('(', ' ')
            line=line.replace(')', ' ')
            line=line.lower()
            line=line.lstrip()

            # if the line starts with use split it and search if we there is an "only" clause and remove everything after that
            if line.startswith("use "):
                use_calls = line.split()

                # remove use open_mp, as it is a library
                for call in exclude_use_list:
                    if call in use_calls:
                        use_calls.remove(call)

                for i, v in enumerate(use_calls):
                    if v=="only":
                        use_calls=use_calls[:i]
                        break

                if "intrinsic" in use_calls:
                    use_calls = []

                # add wahtever is left in the list of use calls to list of module dependencies
                if len(use_calls)>0:
                    for s in use_calls[1:]:
                        file_module_dependence.add(s)

            # we need to look for modules and submodules
            # we will do this separetly as we
            elif line.startswith("module "):
                module_declaration  = line.split()
                if module_declaration[1] not in(module_exclude_keywords):
                    assert len(module_declaration) == 2, "Module declaration line has more than 2 words!!"
                    modules_declared_in_file.add(module_declaration[1])

            # for submodules we have a syntax of the type submodule(module) submodulename
            elif line.startswith("submodule "):
                module_declaration  = line.split()
                assert len(module_declaration) == 3, "SubModule declaration line has more than 3 words!!"
                module_submodule_dependence.add((module_declaration[1],module_declaration[2]))
                submodules_declared_in_file.add(module_declaration[2])

    return file_module_dependence, modules_declared_in_file, module_submodule_dependence, submodules_declared_in_file

def remove_files_from_list(dep_list, files_to_remove):
    for file in files_to_remove:
        if file in dep_list:
            dep_list.pop(file)
    return dep_list

def invert_dictionary(original_dict):
    """
    This method inverts a dictionary, i.e. it returns a new dictionasy whose keys
    are teh values of the orignal dictionary and viceversa. There is no guarantee
    the number of keys of the inverted dictionary is equal to the number of keys
    of the original dictionary:

    input:
        - original_dict: the dictionary to invert. This method assumes the values
                        in the value of the original dictionary is an iterable object

    output:
        - inverted_dict: the inverted dictionary. The for each key the value is a list
    """
    inverted_dict = dict()
    for k, values in original_dict.items():
        for v in values:
            if v in inverted_dict.keys():
                inverted_dict[v].append(k)
            else:
                inverted_dict[v] = [k]
    return inverted_dict

def order_modules_for_compiler(files_to_used_modules_dict, files_to_decleared_module_dict, \
                               files_to_module_submodule_dict, files_to_decleared_submodule_dict):
    """
    This is a recursive routine that given two inputs, i.e. files_to_used_modules_dict and
    files_to_decleared_module_dict (see below), it generates a dependency tree. It does this
    recursively starting from the bottom of the tree, i.e. from source files without dependencies,
    and works its way up until it finds the main file.
    This routine errors out if it cannot find any dependencies for a file in the list.

    inputs:
        - files_to_used_modules_dict: dictionary associating a file to all the modules it uses.
          Keys are a filenames and values are the modules used by each file.
          Values are an iterable object

        - files_to_decleared_module_dict: dictionary associating a file to all the modules decleared in it.
          Keys are a filenames and values are the modules decleared within the file.
          Values are an iterable object

    output:
        - dependency_tree: orderd list of lists. Each list within the list contains
          all the modules that can be compiled after compiling all the modules in the
          previous lists. For example, dependency_tree = [[A.f90],[B.f90, C.f90]]
          means that before compiling B.f90, C.f90 we need to compile A.f90.

    """

    # first we invert to files_to_decleared_module_dict to have a map between modules and
    # the file where they are decleared
    declared_module_to_file_dict = invert_dictionary(files_to_decleared_module_dict)
    declared_submodule_to_file_dict = invert_dictionary(files_to_decleared_submodule_dict)

    # generate a complete file list this will become handy later
    complete_file_list = set()
    for file, _ in files_to_decleared_module_dict.items():
        complete_file_list.add(file)

    # list contatinig the ordered dependency tree
    dependency_tree = list()

    # the set of files that can be compiled. This is updated at every iteration
    compiled_files = set()
    module_files_list = set()
    submodule_files_list = set()
    program_file_list= set()

    # the set of compiled modules and submodules
    compiled_modules_and_submodules = set()

    # we start by compiling all modules. If a file doesn't decelare a module, it either
    # contains a submodule or it is the main program.

    # let's start dividing files into three categories:
    # 1-module files
    # 2-sub module files
    # 3-program files

    #
    for file in complete_file_list:
        if len(files_to_decleared_module_dict[file]) == 0 and \
           len(files_to_decleared_submodule_dict[file]) == 0:
            program_file_list.add(file)
        elif len(files_to_decleared_module_dict[file]) > 0:
            module_files_list.add(file)
        elif len(files_to_decleared_module_dict[file]) == 0 and \
            len(files_to_decleared_submodule_dict[file]) > 0 :
            submodule_files_list.add(file)
        else :
            assert 1<0, "we should not be here, This means a file cannot eb classified"

    assert len(program_file_list) > 0, \
    "The list of program files is empty. This probably means you are declaring \n\
     a module and or a submdoule in the main program. While this is allowed, by \n\
     the fortran standrads it is bad programming practive and this software does NOT \n\
     allow it. Please go back and change your code."

    # remove files not included in module_files_list from files_to_used_modules_dict
    entry_to_remove = set()
    for file, _ in files_to_used_modules_dict.items():
        if file not in module_files_list:
            entry_to_remove.add(file)
    files_to_used_modules_dict = remove_files_from_list(files_to_used_modules_dict, entry_to_remove)

    # recursively search for module dependencies. At each iteration we identify all the
    # compilable files, i.e. the files depneding only from the files in the compiled_files list .
    while len(files_to_used_modules_dict) > 0:

        # the compilable files at this iteration
        compilable_files = list()

        # loop over all the remaining files and find the ones we can compile
        for file, modules in files_to_used_modules_dict.items():
            # list of files in which the modules used by the current file are decleared
            ancestor_files = set()

            # find all the files where the depending modules have been decleared
            for m in modules:
                module_file = declared_module_to_file_dict[m][0]
                if module_file != file:
                    ancestor_files.add(module_file)

            # check if all depending files are already presenet in the compiled_files list
            if ancestor_files.issubset(compiled_files):
                # if so add the current file to the list of compilable files
                compilable_files.append(file)

        # check that at the end of each iteration we identified at list one compilable file
        # if this is not the case teh dependency are broken and the program shall scream an die
        if len(compilable_files) == 0:
            print("Can't find any file that can be compiled!!!")
            print("Remaining files and assocaite modules are: \n")
            for file, modules in files_to_used_modules_dict.items():
                 print("file:", file, "modules:", modules)
            print("\n")
        assert len(compilable_files) > 0, "Can't find any file that can be compiled!!!"

        # add the compilable files list to the full dependency tree
        dependency_tree.append(compilable_files)

        # add the compilable files to the compiled file
        for file in compilable_files:
            compiled_files.add(file)
            # at each pass we will also keep track of wich modules and submodules have been compiled.
            # this will come handy later when we need to compile submodules.
            for decleared_module in files_to_decleared_module_dict[file]:
                compiled_modules_and_submodules.add(decleared_module)

            # check if the compiled file also has submodules
            if file in files_to_decleared_submodule_dict :
                for decleared_submodule in files_to_decleared_submodule_dict[file]:
                    compiled_modules_and_submodules.add(decleared_submodule)
                # if the file also contains submodules, we will remove it from the
                # files_to_decleared_submodule_dictnad from the files_to_module_submodule_dict
                files_to_decleared_submodule_dict.pop(file)
                files_to_module_submodule_dict.pop(file)

        # remove the compilable files from the files_to_used_modules_dict (i.e. from the dictionary we are searching in)
        files_to_used_modules_dict = remove_files_from_list(files_to_used_modules_dict, compilable_files)

    # we are done compiling modules. now we can search the dependecy between modules/submodules and submodules
    # before starting let's remove from files_to_module_submodule_dict all the
    # files that are already been compiled (submodule in those files would have been already compiled)
    # and all the files that are not in the submodule_files_list (e.g. program files)

    entry_to_remove = set()
    for file, _ in files_to_module_submodule_dict.items():
        if (file in compiled_files) or (file not in submodule_files_list):
            entry_to_remove.add(file)
    files_to_module_submodule_dict = remove_files_from_list(files_to_module_submodule_dict, entry_to_remove)

    while len(files_to_module_submodule_dict) > 0:
        compilable_files = list()
        for file, mods_to_subs in files_to_module_submodule_dict.items():
            # find in which file ancestor module and submodules are declared
            ancestor_files =set()
            for mod_to_sub in mods_to_subs:
                ancestor_mod_name = mod_to_sub[0]
                # the ancestor can eihter be a module or a submodule hence we check
                if ancestor_mod_name in declared_submodule_to_file_dict:
                    ancestor_files.add(declared_submodule_to_file_dict[ancestor_mod_name][0])
                elif ancestor_mod_name in declared_module_to_file_dict:
                    ancestor_files.add(declared_module_to_file_dict[ancestor_mod_name][0])
                else:
                    assert 1 < 0, "cannot find the file where the module/submodule named {} is decleared".format(ancestor_mod_name)

            # if the current file is in the ancestor list it means multiple submodules
            # are declared in cascade in the current file, hence we can remove it from the list
            if file in ancestor_files:
                ancestor_files.pop(file)

            # now compare teh list of ancestors with the list of compiled files
            if ancestor_files.issubset(compiled_files):
                # if so add the current file to the list of compilable files
                compilable_files.append(file)

        # check that at the end of each iteration we identified at list one compilable file
        # if this is not the case teh dependency are broken and the program shall scream an die
        if len(compilable_files) == 0:
            print("Can't find any file that can be compiled!!!")
            print("Remaining files and assocaite modules are: \n")
            for file, modules in files_to_used_modules_dict.items():
                 print("file:", file, "modules:", modules)
            print("\n")
        assert len(compilable_files) > 0, "Can't find any file that can be compiled!!!"

        # add the compilable files list to the full dependency tree
        dependency_tree.append(compilable_files)

        # add the compilable files to the compiled file
        for file in compilable_files:
            compiled_files.add(file)
            # at each pass we will also keep track of wich modules and submodules have been compiled.
            # this will come handy later when we need to compile submodules.
            for decleared_submodule in files_to_decleared_submodule_dict[file]:
                compiled_modules_and_submodules.add(decleared_submodule)

        # remove the compilable files from the files_to_used_modules_dict (i.e. from the dictionary we are searching in)
        files_to_module_submodule_dict = remove_files_from_list(files_to_module_submodule_dict, compilable_files)

    # we are done with submodules, at this point we should be able to compile the main files safely
    dependency_tree.append(program_file_list)

    return dependency_tree

def create_compilation_string(dependency_tree, files_path_dict):
    """
    This routine creates the string with filenames that can be used by the compiler to compile the
    fortran90 program.
    inputs:
        - dependency_tree: orderd list of lists. Each list within the list contains
          all the modules that can be compiled after compiling all the modules in the
          previous lists. For example, dependency_tree = [[A.f90],[B.f90, C.f90]]
          means that before compiling B.f90, C.f90 we need to compile A.f90.

        - files_path_dict: a dictionary assocaiting to each file its path

    outputs:
        - s: a string of filenames orderd to comply with the depndency tree.
    """

    s=""
    for files in dependency_tree:
        for f in files:
            path = files_path_dict[f]

            if path != ".":
                s+= os.path.join(path,f)  + " "
            else:
                s+= f + " "
    return s

def generate_fortran90_compiler_script(exclude_use_list=list()):
    """
    Generates the ordered filelist to be used by a fortran90Compiler
    inputs:
        exclude_use_list: The list of used modules to skip (i.e. [lib_1, lib2])
    """
    # find all f90 files within the current working directory
    files_dir_dict = find_F90_files()
    # initialize the dictionary associating files with the modules it uses
    files_to_used_modules = dict()
    # initialize the dictionary associating files with the modules decleared in it
    files_to_decleared_modules = dict()

    # initialize the dictionary associating a file to the pairs <module, submodule>
    files_to_module_submodule = dict()
    # initialize the dictionary associating files with the submodules decleared in it
    files_to_decleared_submodules = dict()

    # search for delcereaed modules and used module in each file
    for file, dir in files_dir_dict.items():
        files_to_used_modules[file], files_to_decleared_modules[file], \
        files_to_module_submodule[file], files_to_decleared_submodules[file] = \
        search_for_USE_and_MODULE_in_file(os.path.join(dir,file), exclude_use_list=exclude_use_list)

    # identify dependencies
    ordered_files_compilation_list = \
    order_modules_for_compiler(files_to_used_modules, files_to_decleared_modules, \
                               files_to_module_submodule, files_to_decleared_submodules)

    # create the fortran string
    return create_compilation_string(ordered_files_compilation_list, files_dir_dict)

if __name__ == "__main__":
    import argparse
    description="\n\
     Runs the generate_fortran90_compiler_script which generates a  \n\
     string of filenames ordered to such that dependencies are respected. \n\
     USE statements of library modules are not part of the main program, and should be skipped \n\
     by adding their name to the  --exclude-use optional argument \n\
     For additional help, open a python shell and execute: \n\
     \n\
     \t import fortran_dependency_resolution \n\
     \t help(fortran_dependency_resolution) \n\
     \n\
     "
    # parsing command linea arguments
    parser=argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--exclude-use', help='A space separated strinf of used modules to skip: e.g. "lib_a lib_b"')
    args=parser.parse_args()

    # populate exclude list
    exclude_use_list = list()
    if args.exclude_use:
        exclude_use_list = args.exclude_use.split()

    # generate the
    print(generate_fortran90_compiler_script(exclude_use_list=exclude_use_list))
