function passed = areEqual(value1, value2, tolerance)
    difference = value1 - value2;
    passed = max(difference, [], 'all') < tolerance;
end