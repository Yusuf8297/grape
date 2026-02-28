# -*- coding: utf-8 -*-
"""
Created on Tue May 10 06:53:28 2022

@author: allan
"""

import re
import math
from operator import attrgetter
import numpy as np
import random
import copy

class Individual(object):
    """
    A GE individual.
    """

    def __init__(self, genome, grammar, max_depth, codon_consumption):
        """
        """
        
        self.genome = genome
        if codon_consumption == 'lazy':
            self.phenotype, self.nodes, self.depth, \
            self.used_codons, self.invalid, self.n_wraps, \
            self.structure = mapper_lazy(genome, grammar, max_depth)
        elif codon_consumption == 'eager':
            self.phenotype, self.nodes, self.depth, \
            self.used_codons, self.invalid, self.n_wraps, \
            self.structure = mapper_eager(genome, grammar, max_depth)
        else:
            raise ValueError("Unknown mapper")

class Grammar(object):
    """
    Attributes:
    - non_terminals: list with each non-terminal (NT);
    - start_rule: first non-terminal;
    - production_rules: list with each production rule (PR), which contains in each position:
        - the PR itself as a string
        - 'non-terminal' or 'terminal'
        - the arity (number of NTs in the PR)
        - production choice label
        - True, if it is recursive, and False, otherwise
        - the minimum depth to terminate the mapping of all NTs of this PR
    - n_rules: df
    
    """
    def __init__(self, file_address):
        #Reading the file
        with open(file_address, "r") as text_file:
            bnf_grammar = text_file.read()
        #Getting rid of all the duplicate spaces
        bnf_grammar = re.sub(r"\s+", " ", bnf_grammar)

        #self.non_terminals = ['<' + term + '>' for term in re.findall(r"\<(\w+)\>\s*::=",bnf_grammar)]
        self.non_terminals = ['<' + term + '>' for term in re.findall(r"\<([\(\)\w,-.]+)\>\s*::=",bnf_grammar)]
        self.start_rule = self.non_terminals[0]
        for i in range(len(self.non_terminals)):
            bnf_grammar = bnf_grammar.replace(self.non_terminals[i] + " ::=", "  ::=")
        rules = bnf_grammar.split("::=")
        del rules[0]
        rules = [item.replace('\n',"") for item in rules]
        rules = [item.replace('\t',"") for item in rules]
        
        #list of lists (set of production rules for each non-terminal)
        self.production_rules = [i.split('|') for i in rules]
        for i in range(len(self.production_rules)):
            #Getting rid of all leading and trailing whitespaces
            self.production_rules[i] = [item.strip() for item in self.production_rules[i]]
            for j in range(len(self.production_rules[i])):
                #Include in the list the PR itself, NT or T, arity and the production choice label
                #if re.findall(r"\<(\w+)\>",self.production_rules[i][j]):
                if re.findall(r"\<([\(\)\w,-.]+)\>",self.production_rules[i][j]):                    
                    #arity = len(re.findall(r"\<(\w+)\>",self.production_rules[i][j]))
                    arity = len(re.findall(r"\<([\(\)\w,-.]+)\>",self.production_rules[i][j]))
                    self.production_rules[i][j] = [self.production_rules[i][j] , "non-terminal", arity, j]
                else:
                    self.production_rules[i][j] = [self.production_rules[i][j] , "terminal", 0, j] #arity 0
        #number of production rules for each non-terminal
        self.n_rules = [len(list_) for list_ in self.production_rules]
  
        for i in range(len(self.production_rules)):
            for j in range(len(self.production_rules[i])):
                NTs_to_check_recursiveness = re.findall(r"\<([\(\)\w,-.]+)\>", self.production_rules[i][j][0])
                NTs_to_check_recursiveness = ['<' + item_ + '>' for item_ in NTs_to_check_recursiveness]
                unique_NTs = np.unique(NTs_to_check_recursiveness, return_counts=False) 
                recursive = False
                for NT_to_check in unique_NTs:
                    stack = [self.non_terminals[i]]  
                    if NT_to_check in stack:
                        recursive = True
                        break
                    else:
                        stack.append(NT_to_check)
                        recursive = check_recursiveness(self, NT_to_check, stack)
                        if recursive:
                            break
                        stack.pop()
                self.production_rules[i][j].append(recursive)
      
        #minimum depth from each non-terminal to terminate the mapping of all symbols
        NT_depth_to_terminate = [None]*len(self.non_terminals)
        #minimum depth from each production rule to terminate the mapping of all symbols
        part_PR_depth_to_terminate = list() #min depth for each non-terminal or terminal to terminate
        isolated_non_terminal = list() #None, if the respective position has a terminal
        #Separating the non-terminals within the same production rule
        for i in range(len(self.production_rules)):
            part_PR_depth_to_terminate.append( list() )
            isolated_non_terminal.append( list() )
            for j in range(len(self.production_rules[i])):
                part_PR_depth_to_terminate[i].append( list() )
                isolated_non_terminal[i].append( list() )
                if self.production_rules[i][j][1] == 'terminal':
                    isolated_non_terminal[i][j].append(None)
                    part_PR_depth_to_terminate[i][j] = 1
                    if not NT_depth_to_terminate[i]:
                        NT_depth_to_terminate[i] = 1
                else:
                    for k in range(self.production_rules[i][j][2]): #arity
                        part_PR_depth_to_terminate[i][j].append( list() )
                        #term = re.findall(r"\<(\w+)\>",self.production_rules[i][j][0])[k]
                        term = re.findall(r"\<([\(\)\w,-.]+)\>",self.production_rules[i][j][0])[k]
                        isolated_non_terminal[i][j].append('<' + term + '>')
        continue_ = True
        while continue_:
            #after filling up NT_depth_to_terminate, we need to run the loop one more time to
            #fill up part_PR_depth_to_terminate, so we check in the beginning
            if None not in NT_depth_to_terminate:
                continue_ = False 
            for i in range(len(self.non_terminals)):
                for j in range(len(self.production_rules)):
                    for k in range(len(self.production_rules[j])):
                        for l in range(len(isolated_non_terminal[j][k])):
                            if self.non_terminals[i] == isolated_non_terminal[j][k][l]:
                                if NT_depth_to_terminate[i]:
                                    if not part_PR_depth_to_terminate[j][k][l]:
                                        part_PR_depth_to_terminate[j][k][l] = NT_depth_to_terminate[i] + 1
                                        if [] not in part_PR_depth_to_terminate[j][k]:
                                            if not NT_depth_to_terminate[j]:
                                                NT_depth_to_terminate[j] = part_PR_depth_to_terminate[j][k][l]
        PR_depth_to_terminate = []
        for i in range(len(part_PR_depth_to_terminate)):
            for j in range(len(part_PR_depth_to_terminate[i])):
                #the min depth to terminate a PR is the max depth within the items of that PR
                if type(part_PR_depth_to_terminate[i][j]) == int:
                    depth_ = part_PR_depth_to_terminate[i][j]
                    PR_depth_to_terminate.append(depth_)
                    self.production_rules[i][j].append(depth_)
                else:
                    depth_ = max(part_PR_depth_to_terminate[i][j])
                    PR_depth_to_terminate.append(depth_)
                    self.production_rules[i][j].append(depth_)
        
def check_recursiveness(self, NT, stack):
    idx_NT = self.non_terminals.index(NT)
    for j in range(len(self.production_rules[idx_NT])):
        NTs_to_check_recursiveness = re.findall(r"\<([\(\)\w,-.]+)\>", self.production_rules[idx_NT][j][0])
        NTs_to_check_recursiveness = ['<' + item_ + '>' for item_ in NTs_to_check_recursiveness]
        unique_NTs = np.unique(NTs_to_check_recursiveness, return_counts=False) 
        recursive = False
  #      while unique_NTs.size and not recursive:
        for NT_to_check in unique_NTs:
            if NT_to_check in stack:
                recursive = True
                return recursive
            else:
                stack.append(NT_to_check) #Include the current NT to check it recursively
                recursive = check_recursiveness(self, NT_to_check, stack)
                if recursive:
                    return recursive
                stack.pop() #If the inclusion didn't show recursiveness, remove it before continuing
    return recursive

def selLexicaseFilterCount(individuals, k):
    """
   

    """
    selected_individuals = []
    #valid_individuals = individuals#.copy()#[i for i in individuals if not i.invalid]
    l_samples = np.shape(individuals[0].fitness_each_sample)[0]
    
    inds_fitness_zero = [ind for ind in individuals if ind.fitness.values[0] == 0]
    if len(inds_fitness_zero) > 0:
        for i in range(k):
            selected_individuals.append(random.choice(inds_fitness_zero))
        return selected_individuals
    
    cases = list(range(0,l_samples))
    candidates = individuals
    
    error_vectors = [ind.fitness_each_sample for ind in candidates]

    unique_error_vectors = list(set([tuple(i) for i in error_vectors]))
    unique_error_vectors = [list(i) for i in unique_error_vectors]
    
    candidates_prefiltered_set = []
    for i in range(len(unique_error_vectors)):
        cands = [ind for ind in candidates if ind.fitness_each_sample == unique_error_vectors[i]]
        candidates_prefiltered_set.append(cands) #list of lists, each one with the inds with the same error vectors

    for i in range(k):
        #fill the pool only with candidates with unique error vectors
        pool = []
        for list_ in candidates_prefiltered_set:
            pool.append(random.choice(list_)) 
        random.shuffle(cases)
        count_ = 0
        while len(cases) > 0 and len(pool) > 1:
            count_ += 1
            f = max
            best_val_for_case = f(map(lambda x: x.fitness_each_sample[cases[0]], pool))
            pool = [ind for ind in pool if ind.fitness_each_sample[cases[0]] == best_val_for_case]
            del cases[0]                    

        pool[0].n_cases = count_
        selected_individuals.append(pool[0]) #Select the remaining candidate
        cases = list(range(0,l_samples)) #Recreate the list of cases

    return selected_individuals
        
def mapper(genome, grammar, max_depth):
    
    idx_genome = 0
    phenotype = grammar.start_rule
    next_NT = re.search(r"\<(\w+)\>",phenotype).group()
    n_starting_NTs = len([term for term in re.findall(r"\<(\w+)\>",phenotype)])
    list_depth = [1]*n_starting_NTs #it keeps the depth of each branch
    idx_depth = 0
    nodes = 0
    structure = []
    
    while next_NT and idx_genome < len(genome):
        NT_index = grammar.non_terminals.index(next_NT)
        index_production_chosen = genome[idx_genome] % grammar.n_rules[NT_index]
        structure.append(index_production_chosen)
        phenotype = phenotype.replace(next_NT, grammar.production_rules[NT_index][index_production_chosen][0], 1)
        list_depth[idx_depth] += 1
        if list_depth[idx_depth] > max_depth:
            break
        if grammar.production_rules[NT_index][index_production_chosen][2] == 0: #arity 0 (T)
            idx_depth += 1
            nodes += 1
        elif grammar.production_rules[NT_index][index_production_chosen][2] == 1: #arity 1 (PR with one NT)
            pass        
        else: #it is a PR with more than one NT
            arity = grammar.production_rules[NT_index][index_production_chosen][2]
            if idx_depth == 0:
                list_depth = [list_depth[idx_depth],]*arity + list_depth[idx_depth+1:]
            else:
                list_depth = list_depth[0:idx_depth] + [list_depth[idx_depth],]*arity + list_depth[idx_depth+1:]

        next_ = re.search(r"\<(\w+)\>",phenotype)
        if next_:
            next_NT = next_.group()
        else:
            next_NT = None
        idx_genome += 1
        
    if next_NT:
        invalid = True
        used_codons = 0
    else:
        invalid = False
        used_codons = idx_genome
    
    depth = max(list_depth)
   
    return phenotype, nodes, depth, used_codons, invalid, 0, structure

def mapper_eager(genome, grammar, max_depth):
    """
    Identical to the previous one.
    Solve the names later.
    """    

    idx_genome = 0
    phenotype = grammar.start_rule
    next_NT = re.search(r"\<(\w+)\>",phenotype).group()
    n_starting_NTs = len([term for term in re.findall(r"\<(\w+)\>",phenotype)])
    list_depth = [1]*n_starting_NTs #it keeps the depth of each branch
    idx_depth = 0
    nodes = 0
    structure = []
    
    while next_NT and idx_genome < len(genome):
        NT_index = grammar.non_terminals.index(next_NT)
        index_production_chosen = genome[idx_genome] % grammar.n_rules[NT_index]
        structure.append(index_production_chosen)
        phenotype = phenotype.replace(next_NT, grammar.production_rules[NT_index][index_production_chosen][0], 1)
        list_depth[idx_depth] += 1
        if list_depth[idx_depth] > max_depth:
            break
        if grammar.production_rules[NT_index][index_production_chosen][2] == 0: #arity 0 (T)
            idx_depth += 1
            nodes += 1
        elif grammar.production_rules[NT_index][index_production_chosen][2] == 1: #arity 1 (PR with one NT)
            pass        
        else: #it is a PR with more than one NT
            arity = grammar.production_rules[NT_index][index_production_chosen][2]
            if idx_depth == 0:
                list_depth = [list_depth[idx_depth],]*arity + list_depth[idx_depth+1:]
            else:
                list_depth = list_depth[0:idx_depth] + [list_depth[idx_depth],]*arity + list_depth[idx_depth+1:]

        next_ = re.search(r"\<(\w+)\>",phenotype)
        if next_:
            next_NT = next_.group()
        else:
            next_NT = None
        idx_genome += 1
        
    if next_NT:
        invalid = True
        used_codons = 0
    else:
        invalid = False
        used_codons = idx_genome
    
    depth = max(list_depth)
   
    return phenotype, nodes, depth, used_codons, invalid, 0, structure

def mapper_lazy(genome, grammar, max_depth):
    """
    This mapper is similar to the previous one, but it does not consume codons
    when mapping a production rule with a single option."""
    
    idx_genome = 0
    phenotype = grammar.start_rule
    next_NT = re.search(r"\<(\w+)\>",phenotype).group()
    n_starting_NTs = len([term for term in re.findall(r"\<(\w+)\>",phenotype)])
    list_depth = [1]*n_starting_NTs #it keeps the depth of each branch
    idx_depth = 0
    nodes = 0
    structure = []
    
    while next_NT and idx_genome < len(genome):
        NT_index = grammar.non_terminals.index(next_NT)
        if grammar.n_rules[NT_index] == 1: #there is a single PR for this non-terminal
            index_production_chosen = 0        
        else: #we consume one codon, and add the index to the structure
            index_production_chosen = genome[idx_genome] % grammar.n_rules[NT_index]
            structure.append(index_production_chosen)
            idx_genome += 1
        
        phenotype = phenotype.replace(next_NT, grammar.production_rules[NT_index][index_production_chosen][0], 1)
        list_depth[idx_depth] += 1
        if list_depth[idx_depth] > max_depth:
            break
        if grammar.production_rules[NT_index][index_production_chosen][2] == 0: #arity 0 (T)
            idx_depth += 1
            nodes += 1
        elif grammar.production_rules[NT_index][index_production_chosen][2] == 1: #arity 1 (PR with one NT)
            pass        
        else: #it is a PR with more than one NT
            arity = grammar.production_rules[NT_index][index_production_chosen][2]
            if idx_depth == 0:
                list_depth = [list_depth[idx_depth],]*arity + list_depth[idx_depth+1:]
            else:
                list_depth = list_depth[0:idx_depth] + [list_depth[idx_depth],]*arity + list_depth[idx_depth+1:]

        next_ = re.search(r"\<(\w+)\>",phenotype)
        if next_:
            next_NT = next_.group()
        else:
            next_NT = None
            
        
    if next_NT:
        invalid = True
        used_codons = 0
    else:
        invalid = False
        used_codons = idx_genome
    
    depth = max(list_depth)
   
    return phenotype, nodes, depth, used_codons, invalid, 0, structure
            
def random_initialisation(ind_class, pop_size, bnf_grammar, 
                          min_init_genome_length, max_init_genome_length,
                          max_init_depth, codon_size, codon_consumption,
                          genome_representation):
        """
        
        """
        population = []
        
        for i in range(pop_size):
            genome = []
            init_genome_length = random.randint(min_init_genome_length, max_init_genome_length)
            for j in range(init_genome_length):
                genome.append(random.randint(0, codon_size))
            ind = ind_class(genome, bnf_grammar, max_init_depth, codon_consumption)
            population.append(ind)
            
        if genome_representation == 'list':
            return population
        elif genome_representation == 'numpy':
            for ind in population:
                ind.genome = np.array(ind.genome)
            return population
        else:
            raise ValueError("Unkonwn genome representation")
    
def sensible_initialisation(ind_class, pop_size, bnf_grammar, min_init_depth, 
                            max_init_depth, codon_size, codon_consumption,
                            genome_representation):
        """
        
        """
        #Calculate the number of individuals to be generated with each method
        is_odd = pop_size % 2
        n_grow = int(pop_size/2)
        
        n_sets_grow = max_init_depth - min_init_depth + 1
        set_size = int(n_grow/n_sets_grow)
        remaining = n_grow % n_sets_grow
        
        n_full = n_grow + is_odd + remaining #if pop_size is odd, generate an extra ind with "full"
        
        #TODO check if it is possible to generate inds with max_init_depth
        
        population = []
        #Generate inds using "Grow"
        for i in range(n_sets_grow):
            max_init_depth_ = min_init_depth + i
            for j in range(set_size):
                remainders = [] #it will register the choices
                possible_choices = [] #it will register the respective possible choices
    
                phenotype = bnf_grammar.start_rule
                remaining_NTs = ['<' + term + '>' for term in re.findall(r"\<([\(\)\w,-.]+)\>",phenotype)] #
                depths = [1]*len(remaining_NTs) #it keeps the depth of each branch
                idx_branch = 0 #index of the current branch being grown
                while len(remaining_NTs) != 0:
                    idx_NT = bnf_grammar.non_terminals.index(remaining_NTs[0])
                    total_options = [PR for PR in bnf_grammar.production_rules[idx_NT]]
                    actual_options = [PR for PR in bnf_grammar.production_rules[idx_NT] if PR[5] + depths[idx_branch] <= max_init_depth_]
                    Ch = random.choice(actual_options)
                    phenotype = phenotype.replace(remaining_NTs[0], Ch[0], 1)
                    depths[idx_branch] += 1
                    if codon_consumption == 'eager':
                        remainders.append(Ch[3])
                        possible_choices.append(len(total_options))
                    elif codon_consumption == 'lazy':
                        if len(total_options) > 1:
                            remainders.append(Ch[3])
                            possible_choices.append(len(total_options))
                    
                    if Ch[2] > 1:
                        if idx_branch == 0:
                            depths = [depths[idx_branch],]*Ch[2] + depths[idx_branch+1:]
                        else:
                            depths = depths[0:idx_branch] + [depths[idx_branch],]*Ch[2] + depths[idx_branch+1:]
                    if Ch[1] == 'terminal':
                        idx_branch += 1
                    
                    remaining_NTs = ['<' + term + '>' for term in re.findall(r"\<([\(\)\w,-.]+)\>",phenotype)]
                
                #Generate the genome
                genome = []
                if codon_consumption == 'eager' or codon_consumption == 'lazy':
                    for k in range(len(remainders)):
                        codon = (random.randint(0,1e10) % math.floor(((codon_size + 1) / possible_choices[k])) * possible_choices[k]) + remainders[k]
                        genome.append(codon)
                else:
                    raise ValueError("Unknown mapper")
                    
                #Include a tail with 50% of the genome's size
                size_tail = max(int(0.5*len(genome)), 1) #Tail must have at least one codon. Otherwise, in the lazy approach, when we have the last PR with just a single option, the mapping procces will not terminate.
                for j in range(size_tail):
                    genome.append(random.randint(0,codon_size))
                    
                #Initialise the individual and include in the population
                ind = ind_class(genome, bnf_grammar, max_init_depth_, codon_consumption)
                
                #Check if the individual was mapped correctly
                if remainders != ind.structure or phenotype != ind.phenotype or max(depths) != ind.depth:
                    raise Exception('error in the mapping')
                    
                population.append(ind)    
            
        for i in range(n_full):
            remainders = [] #it will register the choices
            possible_choices = [] #it will register the respective possible choices

            phenotype = bnf_grammar.start_rule
            remaining_NTs = ['<' + term + '>' for term in re.findall(r"\<([\(\)\w,-.]+)\>",phenotype)] #
            depths = [1]*len(remaining_NTs) #it keeps the depth of each branch
            idx_branch = 0 #index of the current branch being grown

            while len(remaining_NTs) != 0:
                idx_NT = bnf_grammar.non_terminals.index(remaining_NTs[0])
                total_options = [PR for PR in bnf_grammar.production_rules[idx_NT]]
                actual_options = [PR for PR in bnf_grammar.production_rules[idx_NT] if PR[5] + depths[idx_branch] <= max_init_depth]
                recursive_options = [PR for PR in actual_options if PR[4]]
                if len(recursive_options) > 0:
                    Ch = random.choice(recursive_options)
                else:
                    Ch = random.choice(actual_options)
                phenotype = phenotype.replace(remaining_NTs[0], Ch[0], 1)
                depths[idx_branch] += 1
                if codon_consumption == 'eager':
                    remainders.append(Ch[3])
                    possible_choices.append(len(total_options))
                elif codon_consumption == 'lazy':
                    if len(total_options) > 1:
                        remainders.append(Ch[3])
                        possible_choices.append(len(total_options))

                if Ch[2] > 1:
                    if idx_branch == 0:
                        depths = [depths[idx_branch],]*Ch[2] + depths[idx_branch+1:]
                    else:
                        depths = depths[0:idx_branch] + [depths[idx_branch],]*Ch[2] + depths[idx_branch+1:]
                if Ch[1] == 'terminal':
                    idx_branch += 1
                
                remaining_NTs = ['<' + term + '>' for term in re.findall(r"\<([\(\)\w,-.]+)\>",phenotype)]
            
            #Generate the genome
            genome = []
            if codon_consumption == 'eager' or codon_consumption == 'lazy':
            	for j in range(len(remainders)):
            		codon = (random.randint(0,1e10) % math.floor(((codon_size + 1) / possible_choices[j])) * possible_choices[j]) + remainders[j]
            		genome.append(codon)
            else:
            	raise ValueError("Unknown mapper")

            #Include a tail with 50% of the genome's size
            if codon_consumption == 'eager' or codon_consumption == 'lazy':
                size_tail = max(int(0.5*len(genome)), 1) #Tail must have at least one codon. Otherwise, in the lazy approach, when we have the last PR with just a single option, the mapping procces will not terminate.
            
            for j in range(size_tail):
                genome.append(random.randint(0,codon_size))
                
            #Initialise the individual and include in the population
            ind = ind_class(genome, bnf_grammar, max_init_depth, codon_consumption)
            
            #Check if the individual was mapped correctly
            if remainders != ind.structure or phenotype != ind.phenotype or max(depths) != ind.depth:
                raise Exception('error in the mapping')
                
            population.append(ind)    
    
        if genome_representation == 'list':
            return population
        elif genome_representation == 'numpy':
            for ind in population:
                ind.genome = np.array(ind.genome)
            return population
        else:
            raise ValueError("Unkonwn genome representation")
class _PITreeNode:
    """
    A node in the derivation tree built during PI Grow initialisation.
    
    This is an internal helper class used only during PI Grow tree 
    construction. It stores enough information to later extract codons
    in left-to-right (mapper) order.
    
    Attributes
    ----------
    nt : str
        Non-terminal symbol at this node (e.g., '<expr>').
    depth : int
        Depth of this node in the derivation tree.
    prod_idx : int or None
        Index of the chosen production rule for this non-terminal.
    n_choices : int
        Total number of production choices available for this non-terminal.
    children : list of _PITreeNode
        Child nodes (one per non-terminal in the chosen production).
    """
    __slots__ = ['nt', 'depth', 'prod_idx', 'n_choices', 'children']
    
    def __init__(self, nt, depth):
        self.nt = nt
        self.depth = depth
        self.prod_idx = None
        self.n_choices = 0
        self.children = []

def _extract_codons_left_to_right(node, codon_consumption):
    """
    Traverse the derivation tree depth-first left-to-right, collecting
    (prod_idx, n_choices) pairs in the order the GRAPE mapper would 
    consume them.
    
    Parameters
    ----------
    node : _PITreeNode
        Root node to start traversal from.
    codon_consumption : str
        'eager' or 'lazy'. In lazy mode, codons are only recorded when
        the non-terminal has more than one production choice.
    
    Returns
    -------
    list of tuple (int, int)
        List of (production_index, number_of_choices) pairs.
    """
    codons = []
    if codon_consumption == 'eager':
        codons.append((node.prod_idx, node.n_choices))
    elif codon_consumption == 'lazy':
        if node.n_choices > 1:
            codons.append((node.prod_idx, node.n_choices))
    for child in node.children:
        codons.extend(_extract_codons_left_to_right(child, codon_consumption))
    return codons

def _grammar_min_derivation_depth(bnf_grammar):
    """
    Return the shortest derivation depth reachable from the start symbol.
    
    This equals the minimum value of production_rules[0][j][5] across all
    productions of the start non-terminal, where index 5 stores the
    minimum depth needed to derive a complete terminal string.
    
    Parameters
    ----------
    bnf_grammar : Grammar
        A parsed BNF grammar object.
    
    Returns
    -------
    int
        The minimum tree depth at which the grammar can produce a 
        complete (terminal-only) sentence from the start symbol.
    """
    start_prods = bnf_grammar.production_rules[0]
    return min(pr[5] for pr in start_prods)


def _pick_forced_production(all_prods, depth_budget):
    """
    Choose a production when depth-forcing is active (analogous to
    building a *full* tree branch).
    
    Strategy:
    1. Prefer recursive productions that fit within the depth budget,
       so that the tree keeps branching deeper.
    2. If no recursive production fits, accept any production that 
       can terminate within the budget.
    3. Last resort: pick the shallowest-terminating production.
    
    Parameters
    ----------
    all_prods : list
        All production rules for the current non-terminal.
    depth_budget : int
        Remaining depth units from the current node to the target 
        depth (i.e. target_depth − node.depth).
    
    Returns
    -------
    production
        A single selected production rule.
    """
    # Recursive choices that can terminate inside the budget
    recursive_ok = [pr for pr in all_prods
                    if pr[4] and pr[5] <= depth_budget]
    if recursive_ok:
        return random.choice(recursive_ok)
    
    # No recursive choice fits — any choice that terminates in time
    any_ok = [pr for pr in all_prods if pr[5] <= depth_budget]
    if any_ok:
        return random.choice(any_ok)
    
    # Budget is exhausted — pick the shallowest termination path
    min_depth = min(pr[5] for pr in all_prods)
    return random.choice([pr for pr in all_prods if pr[5] == min_depth])


def _pick_grow_production(all_prods, depth_budget):
    """
    Choose a production under standard *grow* semantics: any 
    structurally legal production that terminates within the 
    remaining depth budget.
    
    Parameters
    ----------
    all_prods : list
        All production rules for the current non-terminal.
    depth_budget : int
        Remaining depth units from the current node to the target 
        depth.
    
    Returns
    -------
    production
        A single selected production rule.
    """
    fits = [pr for pr in all_prods if pr[5] <= depth_budget]
    if fits:
        return random.choice(fits)
    
    # Nothing fits — shortest termination path as fallback
    min_depth = min(pr[5] for pr in all_prods)
    return random.choice([pr for pr in all_prods if pr[5] == min_depth])


def PI_Grow(ind_class, pop_size, bnf_grammar, min_init_depth,
            max_init_depth, codon_size, codon_consumption,
            genome_representation):
    """
    Position Independent Grow initialisation for Grammatical Evolution.
    
    Creates a population using the PI Grow method (Fagan et al., 2016),
    which differs from Sensible Initialisation (RHH) in three ways:
    
    1. **Position Independent**: Non-terminals are expanded in random 
       order (not left-to-right), producing more diverse tree shapes.
    2. **Aggressive Depth Forcing**: All non-terminal expansions are 
       restricted to recursive (branching) productions until at least 
       one branch of the tree reaches the target depth.  Additionally,
       if the node being expanded is the last remaining gateway to a 
       recursive sub-tree, it is still forced to branch, preventing 
       premature termination before the depth target is met.
    3. **Grow Only**: No Full trees are generated (unlike RHH which 
       uses half Grow, half Full).
    
    The population is ramped across depths starting from one level 
    above the grammar's minimum derivation depth up to max_init_depth.
    This avoids trivially shallow individuals (single terminals) 
    that reduce initial phenotypic diversity.
    
    Implementation Note
    -------------------
    PI Grow expands non-terminals in random order during tree 
    construction, but GRAPE's mapper (both lazy and eager) reads 
    codons left-to-right. To ensure the genome reproduces the 
    intended phenotype when mapped, we use a two-phase approach:
    
    Phase 1: Build an explicit derivation tree with position-
             independent expansion order and aggressive depth 
             forcing.
    Phase 2: Traverse the completed tree left-to-right (depth-first)
             to extract codons in the order the mapper expects.
    
    Parameters
    ----------
    ind_class : class
        Individual class (typically creator.Individual from DEAP).
    pop_size : int
        Population size.
    bnf_grammar : Grammar
        Parsed BNF grammar object.
    min_init_depth : int
        Minimum initialisation tree depth (may be raised internally 
        to avoid depths below the grammar's minimum derivation depth 
        + 1).
    max_init_depth : int
        Maximum initialisation tree depth.
    codon_size : int
        Maximum codon value (e.g., 255).
    codon_consumption : str
        'eager' or 'lazy'. Determines codon consumption strategy.
    genome_representation : str
        'list' or 'numpy'. Determines genome data structure.
    
    Returns
    -------
    list
        Population of initialised individuals.
    
    References
    ----------
    Fagan, D., Fenton, M. and O'Neill, M. (2016). "Exploring Position
    Independent Initialisation in Grammatical Evolution." IEEE Congress
    on Evolutionary Computation (CEC), pp. 5060-5067.
    """
    # ── Determine ramping range ──────────────────────────────────────
    # PonyGE2 ramps from (min_init_depth + 1) to max_init_depth.
    # The +1 skips the shallowest depth level, which can only produce
    # trivially small trees (often single terminals) that degrade
    # initial population diversity.  The grammar floor additionally
    # prevents depths below the grammar's minimum useful derivation.
    grammar_floor = _grammar_min_derivation_depth(bnf_grammar) + 1
    effective_min = max(min_init_depth + 1, grammar_floor)
    effective_min = min(effective_min, max_init_depth)   # safety clamp
    
    n_depth_levels = max_init_depth - effective_min + 1
    per_level      = pop_size // n_depth_levels
    leftover       = pop_size %  n_depth_levels
    
    # Pre-compute: which NTs in the grammar have at least one
    # recursive production?  Used to tag queue entries so the
    # depth-forcing logic can detect when the last recursive 
    # gateway is about to be consumed.
    nt_has_recursive = []
    for nt_idx in range(len(bnf_grammar.non_terminals)):
        nt_has_recursive.append(
            any(pr[4] for pr in bnf_grammar.production_rules[nt_idx]))
    
    population = []
    
    for level_i in range(n_depth_levels):
        target_depth = effective_min + level_i
        n_inds = per_level + (1 if level_i < leftover else 0)
        
        for _ in range(n_inds):
            # ── Phase 1: Build derivation tree (PI + depth forcing) ──
            
            # Identify non-terminals in the start rule
            start_NTs = ['<' + t + '>' for t in 
                         re.findall(r"\<([\(\)\w,-.]+)\>",
                                    bnf_grammar.start_rule)]
            if not start_NTs:
                start_NTs = [bnf_grammar.start_rule]
            
            root_nodes = [_PITreeNode(nt, 1) for nt in start_NTs]
            
            # Pending queue: each entry is (node, is_recursive_NT).
            # The boolean records whether the NT at that node has any
            # recursive production available — this lets us check 
            # whether removing a node from the queue would leave zero
            # recursive gateways, which would risk the tree 
            # terminating before the target depth is reached.
            pending = []
            for rn in root_nodes:
                idx = bnf_grammar.non_terminals.index(rn.nt)
                pending.append((rn, nt_has_recursive[idx]))
            
            deepest_reached = 1
            
            while pending:
                # ── Position-independent: random node selection ──
                pi_idx  = random.randint(0, len(pending) - 1)
                node, _ = pending.pop(pi_idx)
                
                idx_NT   = bnf_grammar.non_terminals.index(node.nt)
                all_prods = bnf_grammar.production_rules[idx_NT]
                depth_budget = target_depth - node.depth
                
                # Does the queue still contain at least one node 
                # whose NT can recurse?
                queue_has_recursive = any(r for _, r in pending)
                
                # ── Decide whether to force depth ────────────────
                # Force recursive (branching) productions when:
                #   (a) the tree has not yet reached target depth, OR
                #   (b) no other recursive gateway remains in the 
                #       queue (prevents premature termination).
                # Once the target depth has been reached AND there 
                # are still recursive gateways queued, switch to 
                # standard grow semantics.
                force = (deepest_reached < target_depth) or \
                        (not queue_has_recursive)
                
                if force:
                    Ch = _pick_forced_production(all_prods, depth_budget)
                else:
                    Ch = _pick_grow_production(all_prods, depth_budget)
                
                # Store the choice in the tree node for later 
                # codon extraction
                node.prod_idx  = Ch[3]
                node.n_choices = len(all_prods)
                
                child_depth = node.depth + 1
                if child_depth > deepest_reached:
                    deepest_reached = child_depth
                
                # Spawn child nodes for every NT in the chosen 
                # production
                if Ch[1] == 'non-terminal':
                    child_NTs = [
                        '<' + t + '>' for t in
                        re.findall(r"\<([\(\)\w,-.]+)\>", Ch[0])]
                    for cnt in child_NTs:
                        child = _PITreeNode(cnt, child_depth)
                        node.children.append(child)
                        c_idx = bnf_grammar.non_terminals.index(cnt)
                        pending.append(
                            (child, nt_has_recursive[c_idx]))
            
            # ── Phase 2: Extract codons in mapper order ──────────
            # Depth-first left-to-right traversal produces codons 
            # in the order GRAPE's mapper will consume them.
            ordered_codons = []
            for root in root_nodes:
                ordered_codons.extend(
                    _extract_codons_left_to_right(
                        root, codon_consumption))
            
            # Build the genome from the ordered codons
            remainders = []
            genome     = []
            for prod_idx, n_choices in ordered_codons:
                remainders.append(prod_idx)
                codon = (random.randint(0, int(1e10))
                         % math.floor(((codon_size + 1) / n_choices))
                         * n_choices) + prod_idx
                genome.append(codon)
            
            # Append a random tail (≥ 1 codon, ~50 % of used 
            # length).  The tail is required for the lazy mapper 
            # to terminate when the final production has only one 
            # option.
            tail_len = max(int(0.5 * len(genome)), 1)
            for _ in range(tail_len):
                genome.append(random.randint(0, codon_size))
            
            # Instantiate the individual and verify correctness
            ind = ind_class(genome, bnf_grammar, target_depth,
                            codon_consumption)
            
            if remainders != ind.structure or ind.invalid:
                raise Exception('PI_Grow error in mapping')
            
            population.append(ind)
    
    if genome_representation == 'list':
        return population
    elif genome_representation == 'numpy':
        for ind in population:
            ind.genome = np.array(ind.genome)
        return population
    else:
        raise ValueError("Unknown genome representation")
            
def crossover_onepoint(parent0, parent1, bnf_grammar, max_depth, codon_consumption, 
                       genome_representation='list', max_genome_length=None):
    """
    
    """
    if parent0.invalid: #used_codons = 0
        possible_crossover_codons0 = len(parent0.genome)
    else:
        possible_crossover_codons0 = min(len(parent0.genome), parent0.used_codons) #in case of wrapping, used_codons can be greater than genome's length
    if parent1.invalid:
        possible_crossover_codons1 = len(parent1.genome)
    else:
        possible_crossover_codons1 = min(len(parent1.genome), parent1.used_codons)

    parent0_genome = parent0.genome.copy()
    parent1_genome = parent1.genome.copy()
    continue_ = True    
    
    while continue_:
        #Set points for crossover within the effective part of the genomes
        point0 = random.randint(1, possible_crossover_codons0)
        point1 = random.randint(1, possible_crossover_codons1)
      
        if genome_representation == 'list':
            #Operate crossover
            new_genome0 = parent0_genome[0:point0] + parent1_genome[point1:]
            new_genome1 = parent1_genome[0:point1] + parent0_genome[point0:]
        else:
            raise ValueError("Only 'list' representation is implemented")
        
        new_ind0 = reMap(parent0, new_genome0, bnf_grammar, max_depth, codon_consumption)
        new_ind1 = reMap(parent1, new_genome1, bnf_grammar, max_depth, codon_consumption)
  
        continue_ = new_ind0.depth > max_depth or new_ind1.depth > max_depth
    
    if max_genome_length:
        if len(new_ind0.genome) > max_genome_length:
            new_ind0.invalid = True
        if len(new_ind1.genome) > max_genome_length:
            new_ind1.invalid = True
        
    del new_ind0.fitness.values, new_ind1.fitness.values
    return new_ind0, new_ind1   

def mutation_int_flip_per_codon(ind, mut_probability, codon_size, bnf_grammar, max_depth, 
                                codon_consumption, max_genome_length=None):
    """

    """
    # Operation mutation within the effective part of the genome
    if ind.invalid: #used_codons = 0
        possible_mutation_codons = len(ind.genome)
    else:
        possible_mutation_codons = min(len(ind.genome), ind.used_codons) #in case of wrapping, used_codons can be greater than genome's length

    continue_ = True
    #genome = ind.genome.copy()
    genome = copy.deepcopy(ind.genome)
    mutated_ = False
    
    while continue_:
        for i in range(possible_mutation_codons):
            if random.random() < mut_probability:
                genome[i] = random.randint(0, codon_size)
                mutated_ = True
    
        new_ind = reMap(ind, genome, bnf_grammar, max_depth, codon_consumption)
        continue_ = new_ind.depth > max_depth
        
    if max_genome_length:
        if len(new_ind.genome) > max_genome_length:
            new_ind.invalid = True

    if mutated_:
        del new_ind.fitness.values
    return new_ind,

def reMap(ind, genome, bnf_grammar, max_tree_depth, codon_consumption):
    #TODO refazer todo o reMap para nao copiar o ind
    #
    #ind = Individual(genome, bnf_grammar, max_tree_depth, codon_consumption)
    ind.genome = genome
    if codon_consumption == 'lazy':
        ind.phenotype, ind.nodes, ind.depth, \
        ind.used_codons, ind.invalid, ind.n_wraps, \
        ind.structure = mapper_lazy(genome, bnf_grammar, max_tree_depth)
    elif codon_consumption == 'eager':
        ind.phenotype, ind.nodes, ind.depth, \
        ind.used_codons, ind.invalid, ind.n_wraps, \
        ind.structure = mapper_eager(genome, bnf_grammar, max_tree_depth)
    else:
        raise ValueError("Unknown mapper")
        
    return ind

def replace_nth(string, substring, new_substring, nth):
    find = string.find(substring)
    i = find != -1
    while find != -1 and i != nth:
        find = string.find(substring, find + 1)
        i += 1
    if i == nth:
        return string[:find] + new_substring + string[find+len(substring):]
    return string

def selTournamentWithoutInvalids(individuals, k, tournsize, fit_attr="fitness"):
    """
    A simple tournament selection, which avoid invalid individuals.
    """
    chosen = []
    valid_individuals = [i for i in individuals if not i.invalid]
    while len(chosen) < k:
        aspirants = random.sample(valid_individuals, tournsize)
        chosen.append(max(aspirants, key=attrgetter(fit_attr)))
    return chosen
