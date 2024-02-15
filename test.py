import doctest
import normal_surfaces
# also import whatever new file we make for the keeling functions we will put into a new file

modules = [normal_surfaces]

def print_results(module, results):
    print(module.__name__ + ':')
    print('   %s failures out of %s tests.' % (results.failed, results.attempted))

def doctest_modules(modules, verbose=False, print_info=True):
    finder = doctest.DocTestFinder()
    failed, attempted = 0, 0
    for module in modules:
        runner = doctest.DocTestRunner(verbose=verbose)
        for test in finder.find(module):
            runner.run(test)
        result = runner.summarize()
        failed += result.failed
        attempted += result.attempted
        if print_info:
            print_results(module, result)

    if print_info:
        print('\nAll doctests:\n   %s failures out of %s tests.' % (failed, attempted))
    return doctest.TestResults(failed, attempted)

if __name__ == '__main__':
    doctest_modules(modules, True)
