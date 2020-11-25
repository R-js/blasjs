/* eslint-disable @typescript-eslint/no-var-requires */
module.exports = {
    preset: 'ts-jest',
    testEnvironment: 'node',
    verbose: true,
    cacheDirectory: '.jest-cache',
    testPathIgnorePatterns: ['es6', 'commonjs'],
    testMatch: ['**/__tests__/**/*.[t]s?(x)', '**/?(*.)+(spec|test).[t]s?(x)'],
    globals: {
        'ts-jest': {
            compiler: 'typescript',
            tsconfig: 'tsconfig-jest.json',
            diagnostics: {
                ignoreCodes: [151001],
            },
        },
    },
};
