const testRegex = [
    'lib/l3/complex/csyrk/__test__/test.ts'
];

const collectCoverageFrom = [
    'lib/l3/**/*.ts',
];

module.exports = {
    automock: false,
    collectCoverage: true,
    maxWorkers: "50%",
    collectCoverageFrom,
    coveragePathIgnorePatterns: ['node_modules', 'test', 'doc.ts'],
    coverageDirectory: 'coverage',
    //coverageProvider: 'babel', //"v8" is still experimental, but use "v8" for walk through debugging
    coverageProvider: 'v8', //"v8" is still experimental, but use "v8" for walk through debugging
    coverageReporters: ['json', 'lcov', 'text', 'clover'],
    preset: 'ts-jest',
    testEnvironment: 'node',
    verbose: true,
    cacheDirectory: '.jest-cache',
    testPathIgnorePatterns: ['/esm/', '/commonjs/', '/types/'],
    //testMatch: ['**/__tests__/**/*.[t]s?(x)', '**/?(*.)+(spec|test).[t]s?(x)'],
    testRegex,
    transform: {
        "\\.test\\.ts$" :["ts-jest", {
            compiler: 'typescript',
            tsconfig: 'tsconfig.json',
            diagnostics: {
                ignoreCodes: [151001],
            },
        }]
    },
    moduleNameMapper: {
        '^@utils/(.*)$': '<rootDir>/lib/utils/$1',
        /*'^@special/(.*)$': '<rootDir>/src/lib/special/$1',
        '^@trig/(.*)$': '<rootDir>/src/lib/trigonometry/$1',
        '^@rng/(.*)$': '<rootDir>/src/lib/rng/$1',
        '^@lib/(.*)$': '<rootDir>/src/lib/$1',
        '^lib/(.*)$': '<rootDir>/src/lib/$1'
        */
    },
    setupFiles: ['<rootDir>/lib/packages/__test__/jest-ext.d.ts'],
    setupFilesAfterEnv: [
        '<rootDir>/lib/packages/__test__/jest-extension.ts',
        '<rootDir>/lib/packages/__test__/mock-of-debug.ts'
    ],
};

